#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 13 08:59:11 2025

@author: LB & JH
"""

from __future__ import annotations
import os
import sys
from ete4 import Tree
import subprocess
import random
import csv
import scripts.OrthogroupSimulation.coalescence as coalescence
from collections import Counter
from collections import defaultdict
import time
from pathlib import Path
import re
import multiprocessing as mp
from typing import Dict, Tuple
import glob
import shutil
import numpy as np
import string
import math
import statistics

THREADS = 8

## multiprocessing housekeeping ##
# placeholder for shared queue
SEQ_Q = None
# cache to map orthogroup to key and seq (so that failed attempts can retry)
SEED_CACHE: Dict[int, Tuple[str, str]] = {}

## Initializer for each worker process
## args is passed  via pool.starmap
def _init_worker(seq_q):
    global SEQ_Q, SEED_CACHE
    SEQ_Q = seq_q
    SEED_CACHE = {}
######################################

### define other required variables for simulations
tag = "XTAGX"

########################
### Parameter and input helpers
########################

## load parameters from parameter file
## Returns a plain dict of key - value.
# every parameter (iqtree_path, sagephy_path, gbc etc) lives on args, rather than global
def LoadParameters(file_path):
    params = {}
    with open(file_path) as f:
        for line in f:
            line = line.strip()
            if not line or '=' not in line:
                continue
            key, value = line.split('=', 1)
            try:
                if '.' in value or 'e' in value.lower():
                    value = float(value)
                else:
                    value = int(value)
            except ValueError:
                if value.lower() in ("true", "yes", "on"):
                    value = True
                elif value.lower() in ("false", "no", "off"):
                    value = False
                else:
                    value = value  # keep as string
            params[key] = value
    return params

## needed if we are going to do shuffling
def FindRootNode(args):
    species_tree_file = str(args.s)
    t1 = ete4.Tree(species_tree_file, parser=1)
    val = (2*len(t1)-2)
    return val

def FindNSpecies(args):
    species_tree_file = str(args.s)
    t1 = Tree(species_tree_file, parser=1)
    return len(t1)

##################################################
# Gene tree simulation and branch handling
##################################################

## check if tree is at root of species tree
## if not, we need to shuffle
def CheckForRoot(outname, args):
    NS = FindNSpecies(args)
    tree_file = os.path.join(str(args.o), 'temporary_files', f"{outname}.unpruned.tree")
    with open(tree_file, "r") as f:
        tree_str = f.read()
    # minimum DISCPT value = 95% of NS (rounded up)
    threshold = math.ceil(NS * 0.90)
    # find all DISCPT values in the tree
    matches = re.findall(r"DISCPT=\((\d+),", tree_str)
    # return 1 if any DISCPT is within threshold–NS range
    for m in matches:
        val = int(m)
        if threshold <= val <= NS:
            return 1
    return 0

# simulate a gene tree using sagephy
def SimulateGeneTree(outname, args):
    # outname is the orthogroup name (e.g. 11)
    # use random to choose values for parameters where user inputs max
    dupes = random.random() * args.max_duplication_rate
    loss = random.random() * args.max_loss_rate
    trans_rate = random.random() * args.max_transfer_rate
    # define output path
    output1 = os.path.join(str(args.o), 'temporary_files', str(outname))
    # command to be run
    cmd = [
    "java", "-jar", str(args.sagephy_path), "GuestTreeGen",
    "-a", "250",
    "-gb", "-gbc", str(args.gbc),
    "-vp", str(tag),
    "-p", str(args.leaf_sampling_probability),
    "-rt", str(args.replacement_prob),
    str(args.s), str(dupes), str(loss), str(trans_rate), str(output1)]
    #print(" ".join(cmd))
    proc = subprocess.run(cmd, check=True, timeout=60, capture_output=True, text=True)
    # SaGePhy can reject its options and still exit normally, writing no tree,
    # so check that the tree is there rather than failing later on a missing file
    if not os.path.isfile(output1 + ".pruned.tree"):
        raise RuntimeError(
            f"SaGePhy GuestTreeGen wrote no gene tree for orthogroup {outname}. "
            "Its output was:\n" + (proc.stdout + proc.stderr)[-1500:]
        )
    # add to logfile for this orthogroup
    logfile = os.path.join(str(args.o), 'temporary_files', f"{outname}.txt")
    with open(logfile, "a") as f:
        f.write(f"orthogroup={outname}\n")
        f.write(f"duplication_rate={dupes}\n")
        f.write(f"loss_rate={loss}\n")
        f.write(f"transfer_rate={trans_rate}\n")
        #f.write(f"indel_size={indel_size}\n")

# read duplication node IDs out of a SaGePhy guest2host file
## same file/format PrepareForOrthologs() reads for SPECIATION nodes, just
## pulling DUPLICATION rows instead. Kept as its own small helper (rather
## than folding into PrepareForOrthologs) so that function's ground-truth
## orthology behaviour stays completely untouched by the ILS work.
def ReadDuplicationIDs(g2h_file):
    dup_ids = set()
    with open(g2h_file, "r") as f:
        for line_num in range(3):
            next(f)
        for line in f:
            if not line.startswith(tag):
                continue
            parts = line.strip().split()
            if len(parts) >= 3 and parts[2].upper().startswith("DUPLICATION"):
                dup_ids.add(parts[1])
    return dup_ids

# ILS is optional and off by default. It is on only when the user has set
# --ils-theta, which then shows up as args.ils_theta.
def ILSEnabled(args):
    return getattr(args, "ils_theta", None) not in (None, "")

# rearrange the simulated gene tree under the multispecies coalescent (MSC),
# introducing incomplete lineage sorting (ILS) at a severity controlled by
# args.ils_theta (tree depth in coalescent units: larger = less ILS).
# Only called when ILSEnabled(args) is true.
#
# IMPORTANT: this only changes the tree used for branch relaxation ->
# relabelling -> sequence simulation. Ground-truth orthology/paralogy
# (PrepareForOrthologs/Orthologs, below) keeps reading the original,
# untouched .pruned.tree + .pruned.guest2host - orthology is a
# locus-tree/duplication-history property and should not move just because
# the gene tree used to simulate sequences becomes topologically discordant.
def ApplyILSCoalescent(outname, args):
    theta = float(args.ils_theta)
    if theta <= 0:
        raise ValueError(
            f"ils_theta must be > 0 (got {theta}). A value of 0 is the "
            "extreme case of maximal ILS (every branch has zero "
            "coalescent-time length), not an off switch."
        )

    tree_file = os.path.join(str(args.o), 'temporary_files', str(outname) + ".pruned.tree")
    g2h_file = os.path.join(str(args.o), 'temporary_files', str(outname) + ".pruned.guest2host")
    out_file = os.path.join(str(args.o), 'temporary_files', str(outname) + ".pruned_ils.tree")

    with open(tree_file, "r") as f:
        raw = f.read().strip()
    clean = re.sub(r"\[.*?\]", "", raw, flags=re.S)
    raw_tree = Tree(clean, parser=1)

    dup_ids = ReadDuplicationIDs(g2h_file)

    # nondeterministic, consistent with how dupes/loss/transfer rates are
    # already drawn (unseeded random.random()) elsewhere in this pipeline
    seed = random.randint(0, 2**32 - 1)
    coalescent_tree, norm_rf = coalescence.main(
        raw_tree, dup_ids, theta=theta, seed=seed, tag=tag
    )

    coalescent_tree.write(parser=1, outfile=out_file)

    # add to logfile for this orthogroup
    logfile = os.path.join(str(args.o), 'temporary_files', f"{outname}.txt")
    with open(logfile, "a") as f:
        f.write(f"ils_theta={theta}\n")
        f.write(f"ils_norm_rf={norm_rf}\n")

# how many times BranchRelaxer is tried (with a new sigma each time) before giving up
MAX_RELAX_ATTEMPTS = 10

# relax branch lengths using sagephy
def RelaxBranchLengths(outname, args):
    # names need to be species_X_Y
    # where X = orthogroup ID (outname)
    # where Y = ID within that species
    # e.g. homo_sapiens_0_2 is the second human gene in that OG
    # outname is the orthogroup name (e.g. 12)
    # use random to choose values for parameters where user inputs max
    '''
    start = random.random() * max_start_rate
    sigma = random.random() * max_sigma2
    '''
    start = args.max_start_rate
    # define input path - the pruned tree from sagephy, or, when ILS is
    # switched on, the coalescent-rearranged version of it
    # (see ApplyILSCoalescent above)
    tree_suffix = ".pruned_ils.tree" if ILSEnabled(args) else ".pruned.tree"
    input2 = os.path.join(str(args.o), 'temporary_files', str(outname)+tree_suffix)
    # define output path - relaxed pruned tree
    output2 = os.path.join(str(args.o), 'temporary_files', str(outname)+".pruned_relaxed.tree")
    # On a very deep tree (for example under very strong ILS) BranchRelaxer
    # can fail to find valid rates for a large sigma. It then exits normally
    # and writes no tree, so we draw sigma again, up to MAX_RELAX_ATTEMPTS times.
    for attempt in range(MAX_RELAX_ATTEMPTS):
        sigma = np.random.lognormal(mean=args.sigma_log_mean, sigma=args.sigma_log_sd) ** 2
        sigma = min(sigma, 100)
        # command to be run
        cmd = [
        "java", "-jar", str(args.sagephy_path), "BranchRelaxer",
        input2, str(args.relax_model), str(start), str(sigma), "--keep-interior-names",
        "-o", output2]
        proc = subprocess.run(cmd, check=True, timeout=60, capture_output=True, text=True)
        if os.path.isfile(output2):
            break
    else:
        raise RuntimeError(
            f"SaGePhy BranchRelaxer wrote no tree for orthogroup {outname} after "
            f"{MAX_RELAX_ATTEMPTS} attempts. Its output was:\n" + (proc.stdout + proc.stderr)[-1500:]
            + ("\nVery small ils_theta values make gene trees very deep, which can cause this."
               if ILSEnabled(args) else "")
        )
    # add to logfile for this orthogroup
    logfile = os.path.join(str(args.o), 'temporary_files', f"{outname}.txt")
    with open(logfile, "a") as f:
        f.write(f"branch_relax_startrate={start}\n")
        f.write(f"branch_relax_sigma={sigma}\n")

def RelabelRelaxedTree(outname, args):
    # use leafmap to relabel
    filename = os.path.join(str(args.o), 'temporary_files', str(outname)+".pruned.leafmap")
    # we will relabel the relaxed tree
    tree_file = os.path.join(str(args.o), 'temporary_files', str(outname)+".pruned_relaxed.tree")
    # we will name it as relabelled
    tree_outfile = os.path.join(str(args.o), 'temporary_files', str(outname)+".pruned_relaxed_relabelled.tree")
    # we will also record map of old-new labels (helps with orthologs and orthogroups)
    mapping_path = os.path.join(str(args.o), 'temporary_files', str(outname)+".relabelmap.txt")
    mappings = []
    # make a dict to store labels
    my_dict = {}
    with open(filename, 'r') as csv_file:
        reader = csv.reader(csv_file, delimiter='\t')
        for col1, col2 in reader:
            my_dict[col1] = col2
    # need a counter to track how many times each species appears
    species_counts = Counter()
    # load tree in ete4 (we used quoted node names in ete3)
    #tree = ete3.Tree(tree_file, format=1, quoted_node_names=True)
    tree = Tree(tree_file, parser=1)
    #print(tree)
    # relabel tree, and get branch lengths
    total_length = 0
    num_branches = 0
    for node in tree.traverse():
        #if not node.is_root():
        if not node.is_root:
            total_length += node.dist
            num_branches += 1
        #if node.is_leaf():
        if node.is_leaf:
            old_name = node.name
            '''
            node.name2 = my_dict.get(node.name)
            species_counts[node.name2] += 1
            new_name = node.name2 + "_" + str(outname) + "_" + str(species_counts[node.name2])
            '''
            NAME2 = my_dict.get(node.name)
            species_counts[NAME2] += 1
            new_name = NAME2 + "_" + str(outname) + "_" + str(species_counts[NAME2])
            ## and continue with old code
            node.name = new_name
            mappings.append((old_name, new_name))
            
    mean_branch_length = total_length / num_branches
    # write relabelled tree
    #tree.write(format=1, outfile = tree_outfile)
    tree.write(parser=1, outfile = tree_outfile)
    # add to logfile for this orthogroup
    # also log number of genes and number of species
    num_genes = sum(species_counts.values())       # total leaves
    num_species = len(species_counts)
    logfile = os.path.join(str(args.o), 'temporary_files', f"{outname}.txt")
    with open(logfile, "a") as f:
        f.write(f"mean_branch_length={mean_branch_length}\n")
        f.write(f"total_branch_length={total_length}\n")
        #f.write(f"tree={tree.write(format=1)}\n")
        f.write(f"num_genes={num_genes}\n")
        f.write(f"num_species={num_species}\n")
    # save map of old to new name
    with open(mapping_path, "w") as f:
        f.write("Old_Name\tNew_Name\n")
        for old_name, new_name in mappings:
            f.write(f"{old_name}\t{new_name}\n")

##################################################
# Sequence selection and PFAM/domain partitioning
##################################################

# load starting fasta
def MakeFastaDict(args):
    # initialize empty dict
    fasta_dict = {}
    # filepath
    filename = args.f
    # Open the FASTA file and read it into a string
    with open(filename, 'r') as file:
        fasta_data = file.read()
        # Split the data into entries (each entry starts with '>')
        entries = fasta_data.split('>')
        # Process each entry
        for entry in entries:
            if entry.strip():  # Check if the entry is not just whitespace
                lines = entry.splitlines()  # Split the entry into lines
                gene_name = lines[0].split()[0]  # The first line is the gene name
                sequence = ''.join(lines[1:])  # Join the remaining lines to form the sequence
                fasta_dict[gene_name] = sequence  # Add to dictionary
    return fasta_dict

## function to shuffle sequences
def ShuffleSequence(seq, k):
    # Break sequence into k-mers
    kmers = [seq[i:i+k] for i in range(0, len(seq), k)]
    # Shuffle the order of the k-mers
    random.shuffle(kmers)
    # Join back into a string
    return "".join(kmers)

# select a sequence
def SelectSequence(outname, args):
    ## pull one key,seq from the shared queue
    ## on retries (e.g. for alignment too gappy), reuse cached seed
    global SEQ_Q, SEED_CACHE
    if outname in SEED_CACHE:
        key, seq = SEED_CACHE[outname]
    else:
        key, seq = SEQ_Q.get()  # unique, pre-sampled in parent
        SEED_CACHE[outname] = (key, seq)
    # perform shuffling
    # do shuffling if root node doesnt appear in tree
    a = CheckForRoot(outname, args)
    shuffled_sequence = "FALSE"
    def random_char(y):
       return ''.join(random.choice(string.ascii_letters) for x in range(y))
    if a == 0:
        seq = ShuffleSequence(seq, 20) # kmer size = 20
        # current bodge so that no domains in these seqs
        shuffled_sequence = "TRUE"
        key = random_char(5) #can give random name if we dont want domain
    # Log selection
    tmpdir = Path(str(args.o)) / "temporary_files"
    tmpdir.mkdir(parents=True, exist_ok=True)
    logfile = tmpdir / f"{outname}.txt"
    with open(logfile, "a") as f:
        f.write(f"starting_gene_id={key}\n")
        f.write(f"starting_sequence={seq}\n")
        f.write(f"seq_length={len(seq)}\n")
        f.write(f"shuffled_sequence={shuffled_sequence}\n")
    # Write root sequence (overwrite on retry)
    rootseqfile = tmpdir / f"{outname}.root.seq"
    with open(rootseqfile, "w") as f:
        f.write(f">{outname}\n{seq}\n")
    return key, seq

# function to extract pfam domains from pfam run
def ExtractPfamDomains(outname, args):
    pfam_file = args.PFAM
    log_path = os.path.join(str(args.o), 'temporary_files', f"{outname}.txt")
    with open(log_path) as fh:
        for line in fh:
            if line.startswith("starting_gene_id"):
                prefix, gene_id = line.split("=", 1)
                gene_id = gene_id.strip()
    results = []
    types = {"Domain", "Motif"}
    with open(pfam_file) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            cols = line.split()
            if cols[0] != gene_id:
                continue
            if types and cols[7] not in types:
                continue
            hit = {
                "seq_id": cols[0],
                "start": int(cols[1]),
                "end": int(cols[2]),
            }
            results.append(hit)
    # merge overlaps
    if results:
        results.sort(key=lambda x: x["start"])
        merged = [results[0]]
        for h in results[1:]:
            last = merged[-1]
            if h["start"] <= last["end"]:  # overlap
                last["end"] = max(last["end"], h["end"])
            else:
                merged.append(h)
        results = merged
    return results

def SelectRandomDomainModel(args):
    csv_file = str(args.d)
    values = []
    with open(csv_file, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            if row.get("alisim_model"):
                values.append(row["alisim_model"])
    return random.choice(values)


# make partition file from PFAM
def MakePartitionPFAM(outname, args):
   logfile = os.path.join(str(args.o), 'temporary_files', f"{outname}.txt")
   out_path = os.path.join(str(args.o), 'temporary_files', str(outname)+".partionfile.txt")
   root_seq_path = os.path.join(str(args.o), 'temporary_files', f"{outname}.root.seq")
   region_info_path = os.path.join(str(args.o),"temporary_files",str(outname) + ".partition_regions.tsv")
   # set variables
   gamma_shape = random.lognormvariate(mu=args.gamma_shape_mean, sigma=args.gamma_shape_sd)
   prop_invar = random.lognormvariate(mu=args.prop_invar_mean, sigma=args.prop_invar_sd)
   prop_invar = max(0, min(0.99, prop_invar))
   model = f"JTT+G4{{{gamma_shape}}}+I{{{prop_invar}}}"
   model2 = SelectRandomDomainModel(args)
   if model2.startswith("JTT"):
        model2 = model2.replace("JTT", "JTT+C20", 1)
   if model2.startswith("LG"):
        model2 = model2.replace("LG", "LG+C20", 1)
   # get length of seq
   seq_chunks = []
   with open(root_seq_path) as fh:
       for line in fh:
           if line.startswith(">"):
               continue
           seq_chunks.append(line.strip())
   L = len("".join(seq_chunks))
   #
   domain_hits = ExtractPfamDomains(outname, args)
   domain_intervals = [(max(1, h["start"]), min(L, h["end"])) for h in domain_hits if 1 <= h["start"] <= L]
   # Build regions covering the sections, and add FAST gaps and SLOW domains in order.
   regions = []
   cur = 1
   for s, e in sorted(domain_intervals, key=lambda t: t[0]):
       if e < cur:
           # this domain is entirely before the current cursor (due to overlap) = skip
            continue
       if s > cur:
           # gap before this domain = FAST
           regions.append(("gap", cur, s-1, 1, model))
        # add the domain itself = SLOW
       regions.append(("domain", s, e, 0.25, model2))
       cur = e + 1
   if cur <= L:
        # trailing gap to end = FAST
       regions.append(("gap", cur, L, 1, model))
    # If no domains, whole sequence is one FAST region
   if not regions:
        regions = [("gap", 1, L, 1, model)]
    # name regions gene_1..gene_N in order
   regions_named = []
   for i, (region_type, s, e, scalar, m) in enumerate(regions, start=1):
       regions_named.append((f"gene_{i}", region_type, s, e, scalar, m))
   # build NEXUS text
   lines = ["#nexus", "begin sets;"]
   for name, region_type, s, e, _, _ in regions_named:
       lines.append(f"    charset {name} = {s}-{e};")
   lines.append("    charpartition mine = ")
   part_lines = [
    f"        {m}:{name}{{{scalar}}}"
    for name, region_type, _, _, scalar, m in regions_named]
   lines.append(",\n".join(part_lines) + ";")
   lines.append("end;")
   text = "\n".join(lines)
   with open(out_path, "w") as fh:
        fh.write(text)
   with open(logfile, "a") as f:
       f.write(f"gamma_shape_model1={gamma_shape}\n")
       f.write(f"proportion_invariant_sites_model1={prop_invar}\n")
       f.write(f"model2={model2}\n")
       if all(m == model for _,_,_,_,m in regions):
           f.write("models_used=1\n")
       elif all(m == model2 for _,_,_,_,m in regions):
           f.write("models_used=2\n")
       else:
            f.write("models_used=both\n")
   with open(region_info_path, "w") as out:
       out.write("region_name\tregion_type\tstart\tend\n")
       for name, region_type, s, e, _, _ in regions_named:
           out.write(f"{name}\t{region_type}\t{s}\t{e}\n")

##################################################
# Alignment simulation and quality checks
##################################################

# simulate an alignment using edge-proportional partitioning
def SimulateAlignmentPartition(outname, args):
    # outname is the orthogroup name (e.g. 11)
    indel_insert = random.random() * args.max_indel_insert
    indel_delete = random.random() * args.max_indel_delete
    tree_outfile = os.path.join(str(args.o), 'temporary_files', str(outname)+".pruned_relaxed_relabelled.tree")
    root_seq_path = os.path.join(str(args.o), 'temporary_files', f"{outname}.root.seq")
    partition_path = os.path.join(str(args.o), 'temporary_files', str(outname)+".partionfile.txt")
    #indel_size = "GEO{2.5},GEO{2.5}"
    indel_size = "POW{1.7/50},POW{2.2/40}"
    #indel_size = "GEO{1},GEO{1}"
    # define output path
    output1 = os.path.join(str(args.o), 'temporary_files', str(outname) + ".alignment")
    # command to be run
    cmd = [
    str(args.iqtree_path), "--quiet", "--alisim", str(output1),
    "--indel", f"{indel_insert},{indel_delete}",
    "--indel-size", f"{indel_size}",
    "-p", str(partition_path),
    "-t", str(tree_outfile),
    "--out-format", "fasta",
    "--root-seq", f"{root_seq_path},{outname}"
    ]
    logfile = os.path.join(str(args.o), 'temporary_files', f"{outname}.txt")
    with open(logfile, "a") as f:
        f.write(f"indel_insert={indel_insert}\n")
        f.write(f"indel_delete={indel_delete}\n")
        f.write(f"indel_size={indel_size[indel_size.find('{')+1 : indel_size.find('}')]}\n")
    #subprocess.run(cmd, check=True, timeout=60)
    for attempt in range(5):
        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=120
            )
            output = result.stdout + result.stderr
            alignment_file = output1 + ".fa"
            if (
                result.returncode == 0
                and "ERROR:" not in output
                and os.path.isfile(alignment_file)
            ):
                return output
            print(f"Attempt {attempt+1} failed with output:\n{output}")
        except subprocess.TimeoutExpired:
            print(f"Attempt {attempt+1} timed out for orthogroup {outname}")
        time.sleep(1)
    print(f"All AliSim retries failed for orthogroup {outname}!")
    return None    

## function that checks alignment (seqs that are gaps)
def CheckAlignmentHealth(outname, args):
    alignment_path = os.path.join(str(args.o), 'temporary_files', str(outname) + ".alignment.fa")
    if not os.path.isfile(alignment_path):
        return False
    gap_chars="-"
    gapset = set(gap_chars)
    with open(alignment_path) as fh:
        seq = ""
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if seq and set(seq) <= gapset:
                    return True
                seq = ""
            else:
                seq += line
        if seq and set(seq) <= gapset:
            return True
    return False

## function that reads number of duplications and losses
def CountDuplicationLossTransfer(outname, args):
    # duplication and transfer has to be read from pruned (that the tree we see)
    # losses has to be read from unpruned (as by definition they are removed from pruned)
    input_file = os.path.join(str(args.o), 'temporary_files', str(outname) + ".pruned.info")
    input_file2 = os.path.join(str(args.o), 'temporary_files', str(outname) + ".unpruned.info")
    with open(input_file, "r") as f:
        for line in f:
            if line.startswith("No. of duplications:"):
                num_duplications = int(line.split(":")[1].strip())
            elif line.startswith("No. of additive transfers:"):
                num_add_transfers = int(line.split(":")[1].strip())
            elif line.startswith("No. of replacing transfers:"):
                num_replace_transfers = int(line.split(":")[1].strip())
    with open(input_file2, "r") as f:
        for line in f:
            if line.startswith("No. of losses:"):
                num_losses = int(line.split(":")[1].strip())
    total_transfers = num_replace_transfers + num_add_transfers

    logfile = os.path.join(str(args.o), 'temporary_files', f"{outname}.txt")
    with open(logfile, "a") as f:
        f.write(f"num_losses={num_losses}\n")
        f.write(f"num_duplications={num_duplications}\n")
        f.write(f"num_additive_transfers={num_add_transfers}\n")
        f.write(f"num_replacing_transfers={num_replace_transfers}\n")
        f.write(f"total_transfers={total_transfers}\n")

##################################################
# Gap shuffling
##################################################

def ReadGapPositionProfile(profile_file):
    """
    Read gap_position_profile_counts.tsv.
    This function randomly selects one empirical profile row, then returns
    its bins sorted from least gap-enriched to most gap-enriched.
    """
    profiles = []
    with open(profile_file) as f:
        header = next(f).strip().split("\t")

        if "total_gaps" not in header:
            raise ValueError(
                f"Expected a 'total_gaps' column in {profile_file}. "
                "This function a total gaps col"
            )

        total_gaps_idx = header.index("total_gaps")

        bin_columns = [
            col
            for col in header
            if col.startswith("bin_")
        ]

        if len(bin_columns) != 100:
            raise ValueError(
                f"Expected 100 bin columns in {profile_file}, "
                f"but found {len(bin_columns)}."
            )

        # Sort bin columns numerically so bin_10 does not come before bin_2.
        bin_columns = sorted(
            bin_columns,
            key=lambda x: int(x.split("_")[1])
        )

        bin_indices = [
            header.index(col)
            for col in bin_columns
        ]

        for line in f:
            if not line.strip():
                continue

            parts = line.strip().split("\t")

            total_gaps = int(parts[total_gaps_idx])

            # Ignore empirical alignments with no gaps.
            if total_gaps <= 0:
                continue

            probabilities = [
                float(parts[i])
                for i in bin_indices
            ]

            profiles.append(probabilities)

    if not profiles:
        return []

    # Randomly select one empirical alignment profile.
    chosen_profile = random.choice(profiles)

    rows = [
        (b, chosen_profile[b])
        for b in range(100)
    ]

    # Rank bins from least gap-enriched to most gap-enriched.
    rows = sorted(rows, key=lambda x: x[1])

    return rows


def ReadAlignment(alignment_file):
    """
    Read FASTA alignment.
    """
    names = []
    seqs = []

    with open(alignment_file) as f:
        name = None
        seq = []

        for line in f:
            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):
                if name is not None:
                    names.append(name)
                    seqs.append("".join(seq))

                name = line[1:]
                seq = []

            else:
                seq.append(line)

        if name is not None:
            names.append(name)
            seqs.append("".join(seq))

    return names, seqs


def WriteAlignment(names, seqs, alignment_file):
    """
    Write FASTA alignment.
    """
    with open(alignment_file, "w") as out:
        for name, seq in zip(names, seqs):
            out.write(f">{name}\n{seq}\n")


def GetColumnGapFractions(seqs):
    """
    Calculate gap fraction for each alignment column.
    """
    if not seqs:
        return []

    aln_len = len(seqs[0])
    nseq = len(seqs)

    gap_fractions = []

    for col in range(aln_len):
        gaps = 0

        for seq in seqs:
            if col < len(seq) and seq[col] == "-":
                gaps += 1

        gap_fractions.append(gaps / nseq)

    return gap_fractions


def GetGapColumnBlocks(seqs):
    """
    Identify contiguous blocks of alignment columns where at least one
    sequence has a gap. These blocks are treated as approximate shared indel events.
    """
    gap_fractions = GetColumnGapFractions(seqs)

    blocks = []
    in_block = False
    start = None

    for i, gf in enumerate(gap_fractions):
        if gf > 0 and not in_block:
            start = i
            in_block = True

        elif gf == 0 and in_block:
            blocks.append((start, i - 1))
            in_block = False

    if in_block:
        blocks.append((start, len(gap_fractions) - 1))

    return blocks, gap_fractions


def GapShuffleAlignment(outname, args):
    """
    Shuffle gap-containing column blocks across the whole alignment.
    Logic:
        - find contiguous column blocks where at least one sequence has a gap
        - treat each block as a putative indel block
        - rank blocks by how gappy they are
        - rank empirical alignment positions by how gap-enriched they are
        - reorder the gap blocks so more gappy blocks are moved toward
          empirically gap-rich parts of alignments
        - keep non-gap columns in their original relative order
    """
    alignment_file = os.path.join(
        str(args.o),
        "temporary_files",
        f"{outname}.alignment.fa"
    )

    if not os.path.isfile(alignment_file):
        return

    names, seqs = ReadAlignment(alignment_file)

    if not seqs:
        return

    aln_len = len(seqs[0])

    if aln_len == 0:
        return

    # Sanity check: all sequences should have the same aligned length
    if any(len(seq) != aln_len for seq in seqs):
        raise ValueError(f"Alignment has unequal sequence lengths: {alignment_file}")

    profile_rows = ReadGapPositionProfile(args.gap_profile)

    if not profile_rows:
        return

    blocks, gap_fractions = GetGapColumnBlocks(seqs)

    if not blocks:
        logfile = os.path.join(str(args.o), "temporary_files", f"{outname}.txt")
        with open(logfile, "a") as f:
            f.write("gap_shuffled=FALSE_no_gap_blocks\n")
        return

    # Store information for each gap block
    block_info = []

    for start, end in blocks:
        mean_gap = statistics.mean(gap_fractions[start:end + 1])
        block_len = end - start + 1

        block_info.append({
            "start": start,
            "end": end,
            "length": block_len,
            "mean_gap": mean_gap,
        })

    # Rank simulated gap blocks from least gappy to most gappy
    block_info = sorted(block_info, key=lambda x: x["mean_gap"])

    # Empirical bins ranked from least gap-enriched to most gap-enriched
    empirical_bins_ranked = [
        row[0]
        for row in profile_rows
    ]

    # Assign least-gappy blocks to least-gappy empirical bins,
    # and most-gappy blocks to most-gappy empirical bins.
    assigned = []
    n_blocks = len(block_info)

    for rank, block in enumerate(block_info):
        bin_index = int(rank * len(empirical_bins_ranked) / n_blocks)

        if bin_index >= len(empirical_bins_ranked):
            bin_index = len(empirical_bins_ranked) - 1

        target_bin = empirical_bins_ranked[bin_index]

        # random.random() breaks ties within the same target bin
        assigned.append((target_bin, random.random(), block))

    # Sort by empirical target position
    assigned.sort()

    shuffled_blocks = [
        block
        for _, _, block in assigned
    ]

    # Mark all original gap-block columns
    gap_block_cols = set()

    for start, end in blocks:
        for col in range(start, end + 1):
            gap_block_cols.add(col)

    # Non-gap columns are all columns that are not part of a gap block.
    # They will remain in their original relative order.
    non_gap_cols = [
        col
        for col in range(aln_len)
        if col not in gap_block_cols
    ]

    # Rebuild the column order:
    # - walk through the original gap-block slots
    # - copy intervening non-gap columns in original order
    # - replace each original gap block with one shuffled gap block
    new_order = []
    non_gap_idx = 0
    block_idx = 0

    original_blocks_sorted = sorted(blocks)

    for original_start, original_end in original_blocks_sorted:

        # Add non-gap columns before this original gap-block slot
        while (
            non_gap_idx < len(non_gap_cols)
            and non_gap_cols[non_gap_idx] < original_start
        ):
            new_order.append(non_gap_cols[non_gap_idx])
            non_gap_idx += 1

        # Add the next shuffled gap block
        block = shuffled_blocks[block_idx]
        new_order.extend(range(block["start"], block["end"] + 1))
        block_idx += 1

    # Add remaining non-gap columns after the final original gap block
    while non_gap_idx < len(non_gap_cols):
        new_order.append(non_gap_cols[non_gap_idx])
        non_gap_idx += 1

    if len(new_order) != aln_len:
        raise ValueError(
            f"Gap shuffling changed alignment length for {outname}: "
            f"{len(new_order)} vs {aln_len}"
        )

    if len(set(new_order)) != aln_len:
        raise ValueError(
            f"Gap shuffling duplicated or lost columns for {outname}"
        )

    shuffled_seqs = []

    for seq in seqs:
        shuffled_seq = "".join(seq[i] for i in new_order)
        shuffled_seqs.append(shuffled_seq)

    WriteAlignment(names, shuffled_seqs, alignment_file)

    logfile = os.path.join(str(args.o), "temporary_files", f"{outname}.txt")

    with open(logfile, "a") as f:
        f.write("gap_shuffled=TRUE_whole_alignment_gap_blocks\n")
        f.write(f"gap_blocks={len(blocks)}\n")
        f.write(f"gap_block_columns={len(gap_block_cols)}\n")
        f.write(f"non_gap_columns={len(non_gap_cols)}\n")
        f.write(f"alignment_length={aln_len}\n")

##################################################
# Ortholog inference and parsing
##################################################

### Orthologs
## we do this in three parts
## first, we PrepareForOrthologs [clean the tree and prepare a dict of interal nodes]
## then we do our main Ortholog function to call orthologs and save as pairs ###
## then we OrthologParsing to add species names and format the file for that og###

def PrepareForOrthologs(outname, args):
    # prepare tree
    tree_file = os.path.join(str(args.o), "temporary_files", f"{outname}.pruned.tree")
    with open(tree_file, "r") as f:
        raw = f.read().strip()
    clean = re.sub(r"\[.*?\]", "", raw, flags=re.S)
    tree = Tree(clean, parser=1)
    ## guest2host file
    g2h_file = os.path.join(str(args.o), "temporary_files", f"{outname}.pruned.guest2host")
    spec_nodes = set()
    with open(g2h_file, "r") as f:
        for line_num in range(3):
            next(f)
        for line in f:
            if not line.startswith(tag):
                continue
            parts = line.strip().split()
            # we only care about speciation nodes
            if len(parts) >= 3 and parts[2].upper().startswith("SPECIATION"):
                spec_nodes.add(parts[0])
    name_index = {n.name: n for n in tree.traverse() if n.name}
    return tree, spec_nodes, name_index

## function that extracts orthologs
def Orthologs(outname, args):
## to allow multiprocessing, save one file per OG
    ## load the tree and node dict from PrepareForOrthologs
    tree, spec_nodes, name_index = PrepareForOrthologs(outname, args)
    out_path = os.path.join(str(args.o), "temporary_files", f"{outname}.orthologpairs.txt")
    # we now have a tree, and a dict of speciation nodes
    # orthologs are pairs of genes seperated by a speciation node
    # we record pairs we have seen
    seen = set()
    # we prepare to store pairs
    rows = []
    # we loop through speciation nodes
    for sname in spec_nodes:
        # we use the name index from PrepareForOrthologs for the tree node
        node = name_index[sname]
        # we get the children
        children = node.get_children()
        # if its a polytomy - we panic
        if len(children) != 2:
            print("UH OH - ITS A POLYTOMY!")
            quit()
        # we get the children
        left, right = children
        #left_leaves = [l.name for l in left.iter_leaves()]
        #right_leaves = [l.name for l in right.iter_leaves()]
        left_leaves = [l.name for l in left.leaves()]
        right_leaves = [l.name for l in right.leaves()]
        # all pairs are orthologs
        # we use 'seen' to not 'double-count'
        for a in left_leaves:
            for b in right_leaves:
                pair = tuple(sorted((a, b)))
                if pair in seen:
                        continue
                seen.add(pair)
                rows.append(pair)
    # we save a file
    with open(out_path, "w") as w:
        w.write("Gene_1\tGene_2\n")
        for g1, g2 in rows:
            w.write(f"{g1}\t{g2}\n")

# function that extracts species from gene
# useful helper for orthologs and orthogroups
def SpeciesFromGene(gene):
        # genus_species_og_index
        # use rsplit to get genus_species
        # need to adapt if species names have more _
        parts = gene.rsplit("_", 2)
        return parts[0]

## function that extract species names and relabelled gene names of orthologs
def OrthologParsing(outname, args):
## we want headers of Orthogroup	Species_1	Species_2	Gene_1	Gene_2
## relabel map helps us relabel the genes
## orthlogpairs has the orthologs
    # define paths to labelmap, orthologpairs, and output
    map_path = os.path.join(str(args.o), 'temporary_files', str(outname)+".relabelmap.txt")
    pairs_path = os.path.join(str(args.o), 'temporary_files', str(outname)+".orthologpairs.txt")
    out_path = os.path.join(str(args.o), 'temporary_files', str(outname)+".orthologs.txt")
    # Load Old-New mapping
    old2new = {}
    with open(map_path, "r") as f:
        next(f, None)  # skip header
        for line in f:
            line = line.strip()
            if not line:
                continue
            cols = line.split("\t")
            old2new[cols[0]] = cols[1]
    # write final table
    # read from pairs. then parse. then write to outpath
    with open(out_path, "w") as out, open(pairs_path, "r") as pf:
        out.write("Orthogroup\tSpecies_1\tSpecies_2\tGene_1\tGene_2\n")
        next(pf, None)  # skip header
        for line in pf:
            line = line.strip()
            if not line:
                continue
            cols = line.split("\t")
            g1_old, g2_old = cols[0], cols[1]
            # genes come from getting the relabled gene name
            g1_new = old2new.get(g1_old)
            g2_new = old2new.get(g2_old)
            # species come from splitting the gene name
            s1 = SpeciesFromGene(g1_new)
            s2 = SpeciesFromGene(g2_new)
            out.write(f"{outname}\t{s1}\t{s2}\t{g1_new}\t{g2_new}\n")


##################################################
# One-orthogroup and multiprocessing pipeline
##################################################

# below is the old RunOrthogroup that didnt use multiproc
## main function to simulate an orthogroup
def RunOrthogroup(outname, args):
    SimulateGeneTree(outname, args)
    if ILSEnabled(args):
        ApplyILSCoalescent(outname, args)
    RelaxBranchLengths(outname, args)
    RelabelRelaxedTree(outname, args)
    SelectSequence(outname, args)
    MakePartitionPFAM(outname, args)
    alisim_result = SimulateAlignmentPartition(outname, args)
    if alisim_result is None:
        print(f"[{outname}] AliSim failed. Skipping this orthogroup.")
        return
    GapShuffleAlignment(outname, args)
    unhealthy = CheckAlignmentHealth(outname, args)
    ## if check alignment health = true, we need to start again
    ## for that og. remove files with the outname. prefix, and print a message
    if unhealthy:
        print(f"[{outname}] alignment has a sequence with all gaps. Retrying…")
        base_outdir = Path(args.o)
        temp_file = base_outdir / "temporary_files" / f"{outname}.txt"
        if temp_file.exists():
            temp_file.unlink()
        return RunOrthogroup(outname, args)

    if not unhealthy:
        CountDuplicationLossTransfer(outname, args)
        Orthologs(outname, args)
        OrthologParsing(outname, args)
    #print("done an orthogroup")

def RunOrthogroupMultiProc(n, outdir, threads, args):
    # Build starting genome once
    # sample N unique starting genes
    # make a shared queue
    # run the RunOrthogroup function
    # Build once and sample without replacement
    fasta_dict = MakeFastaDict(args)  # {gene_id: sequence}
    if n > len(fasta_dict):
        raise ValueError(f"Need {n} sequences but only have {len(fasta_dict)}")
    samples = random.sample(list(fasta_dict.items()), n)  # [(key, seq), ...]

    # Shared queue (spawn context is safest across OSes)
    ctx = mp.get_context("spawn")
    manager = ctx.Manager()
    seq_q = manager.Queue()
    for kv in samples:
        seq_q.put(kv)

    # Fan out: each worker pulls exactly one seed (cached on retry)
    nproc = threads
    with ctx.Pool(processes=nproc, initializer=_init_worker, initargs=(seq_q,)) as pool:
        pool.starmap(RunOrthogroup, [(i, args) for i in range(n)])


##################################################
# Final output collation
##################################################

### function that makes global ortholog output from individual orthogroup file
##output file is Simulated_orthologs.txt
## headers are Orthogroup	Species_1	Species_2	Gene_1	Gene_2
def ConcatOrthologs(args):
    outdir = args.o
    ## find all ortholog files
    files = sorted(glob.glob(os.path.join(outdir, "temporary_files", "*orthologs.txt"), recursive=True))
    # define path
    out_path = os.path.join(outdir, "Simulated_orthologs.txt")
    # cat files, keeping only header from the first
    wrote_header = False
    with open(out_path, "w") as out:
        for fp in files:
            with open(fp, "r") as f:
                for i, line in enumerate(f):
                    if i == 0:
                        if not wrote_header:
                            out.write(line)   # write header from the first file only
                            wrote_header = True
                        # skip header for subsequent files
                    else:
                        out.write(line)

## function that saves alignments

def CopyAlignments(args):
    outdir = args.o
    files = glob.glob(os.path.join(outdir, "temporary_files", "*alignment.fa"), recursive=True)
    dest = os.path.join(outdir, "alignment_files")
    for src in sorted(files):
        shutil.copy2(src, dest)

## function that saves trees
def CopyTrees(args):
    outdir = args.o
    suffix = ".pruned_relaxed_relabelled.tree"
    files = glob.glob(os.path.join(outdir, "temporary_files", f"*{suffix}"))
    dest = os.path.join(outdir, "tree_files")
    for src in sorted(files):
        fname = os.path.basename(src)
        outname = fname[:-len(suffix)]
        shutil.copy2(src, os.path.join(dest, f"{outname}.tre"))

## function that saves proteomes
def BuildProteomes(args):
    outdir = args.o
    infiles = glob.glob(os.path.join(outdir, "temporary_files", "*alignment.fa"))
    dest = os.path.join(outdir, "proteome_files")
    proteome = defaultdict(list)
    for fp in sorted(infiles):
        # open the alignment. get the genes. find the species. add to the right proteome
        with open(fp) as f:
            # start with gene = none, and a blank 'chunk' to store seq
            # we do this as there is nothing to 'flush' for the first record
            gene = None
            chunks = []
            for line in f:
                if line.startswith(">"):
                    # if starts with >, its the gene ID
                    if gene is not None:
                        # get rid of gaps in alignment
                        seq = "".join(chunks).replace("-", "").strip()
                        # use our function to get species from gene name
                        sp = SpeciesFromGene(gene)
                        # add to proteome
                        proteome[sp].append((gene, seq))
                    gene = line[1:].strip()
                    chunks = []
                    # else, its the sequence
                else:
                    chunks.append(line.strip())
            # flush last record
            if gene is not None:
                seq = "".join(chunks).replace("-", "").strip()
                sp = SpeciesFromGene(gene)
                proteome[sp].append((gene, seq))
    for sp in sorted(proteome):
        outpath = os.path.join(dest, f"{sp}.fa")
        with open(outpath, "w") as w:
            for gene, seq in proteome[sp]:
                w.write(f">{gene}\n{seq}\n")


## function that saves orthogroups
def SaveOrthogroups(args):
    outdir = args.o
    infiles = sorted(glob.glob(os.path.join(outdir, "alignment_files", "*.fa")))
    ## need to dictionary to store which species each gene in an OG belongs to
    og_map = defaultdict(lambda: defaultdict(list))
    all_species = set()
    for fp in infiles:
        # get og ID
        og = os.path.basename(fp).split(".", 1)[0]
        with open(fp) as f:
            for line in f:
                if line.startswith(">"):
                    gene = line[1:].strip()
                    # use our function to get species
                    sp = SpeciesFromGene(gene)
                    # add it to the dict
                    all_species.add(sp)
                    og_map[og][sp].append(gene)
    # save it here
    out_path = os.path.join(outdir, "Simulated_Orthogroups.txt")
    species_cols = sorted(all_species)
    with open(out_path, "w") as w:
        w.write("Orthogroup\t" + "\t".join(species_cols) + "\n")
        # try and sort numerically, not alphabetically
        def og_sort_key(x):
            return (0, int(x))
        for og in sorted(og_map.keys(), key=og_sort_key):
            cells = [",".join(og_map[og].get(sp, [])) for sp in species_cols]
            w.write(og + "\t" + "\t".join(cells) + "\n")

## function that collates stats from logfile to txt file / csv
def BuildOrthogroupStats(args):
    outdir = args.o
    #logs = sorted(glob.glob(os.path.join(outdir, "temporary_files", "*.txt")))
    # I should have named the logfile for each og something else
    # because now we have to horrible regex to find the exact file
    logs_dir = os.path.join(outdir, "temporary_files")
    logs = sorted(
    (os.path.join(logs_dir, fn) for fn in os.listdir(logs_dir)
     if re.fullmatch(r"\d+\.txt", fn)),
    key=lambda p: int(os.path.basename(p)[:-4])
)
    # define the order of variables
    VAR_ORDER = [
        "orthogroup",
        "starting_gene_id",
        "seq_length",
        "shuffled_sequence",
        "num_species",
        "num_genes",
        "duplication_rate",
        "num_duplications",
        "loss_rate",
        "num_losses",
        "transfer_rate",
        "num_replacing_transfers",
        "num_additive_transfers",
        "total_transfers",
        "proportion_invariant_sites_model1",
        "gamma_shape_model1",
        "model2",
        "models_used",
        "indel_insert",
        "indel_delete",
        "indel_size",
        "mean_branch_length",
        "total_branch_length",
        "starting_sequence",
        "branch_relax_startrate",
        "branch_relax_sigma",
        "ils_theta",
        "ils_norm_rf",
    ]

    # Parse all logs into rows
    rows = []
    for fp in logs:
        og = os.path.basename(fp).split(".", 1)[0]
        data = {}
        with open(fp) as f:
            for line in f:
                line = line.strip()
                # var name is before = and value is after
                k, v = line.split("=", 1)
                data[k.strip()] = v.strip()
        rows.append((og, data))
    # function to try and maintain numeric order
    def _key(item):
        og = item[0]
        return (0, int(og))
    rows.sort(key=_key)
    # Write output
    out_path = os.path.join(outdir, "Simulated_Orthogroup_statistics.txt")
    with open(out_path, "w") as w:
        w.write("\t".join(VAR_ORDER) + "\n")
        for _, data in rows:
            w.write("\t".join(data.get(col, "") for col in VAR_ORDER) + "\n")

## function that writes run_log
def WriteRunLog(args):
    outdir = args.o
    log_path = os.path.join(outdir, "Run_log.txt")
    # Collect user input flags
    lines = []
    lines.append(f"# Run started: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    for flag in ["s", "f", "n", "p", "o"]:
        lines.append(f"-{flag} = {getattr(args, flag)}")

    # Write header + flags
    with open(log_path, "w") as w:
        w.write("Orthogroup Simulation Run Log\n")
        w.write("\n# User input flags\n")
        w.write("\n".join(lines) + "\n")

        w.write(f"\n# Parameter file: {os.path.abspath(args.p)}\n")
        with open(args.p, "r") as pf:
            w.write(pf.read())
