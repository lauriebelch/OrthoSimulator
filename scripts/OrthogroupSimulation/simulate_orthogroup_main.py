#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# simulate_orthogroup_main.py

#   1. read_config()
#   2. configure_legacy_args()
#   3. create_subfolders()
#   5. og.RunOrthogroupMultiProc()
#   6. build the simulation output files

import os
import time
import multiprocessing as mp
from pathlib import Path
from types import SimpleNamespace

import scripts.OrthogroupSimulation.orthogroup_simulation_utils as og
from scripts.PFAM import pfam_archive

def read_config(config_path):
    """
    Read the OrthoSimulator-generated output/config.txt
    """
    config = {}
    with open(config_path) as handle:
        for line_number, line in enumerate(handle, start=1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if "=" not in line:
                raise ValueError(
                    f"Bad config line {line_number} in {config_path}: {line}"
                )
            key, value = line.split("=", 1)
            config[key.strip()] = value.strip()
    return config

def require_config(config, key):
    """
    Look for required things in config
    """
    value = config.get(key)
    if value is None or str(value).strip() == "":
        raise ValueError(f"Required value missing from OrthoSimulator config: {key}")
    return value


def resolve_path(path_value, project_root):
    """
    Turn a paths from config.txt into absolute paths
    """
    path = Path(path_value)
    if path.is_absolute():
        return str(path)
    return str((Path(project_root) / path).resolve())


def first_existing_file(candidates, description):
    """
    find the right PFAM file
    """
    for candidate in candidates:
        if candidate and os.path.isfile(candidate):
            return candidate
    tried = "\n".join(f"  {candidate}" for candidate in candidates if candidate)
    raise ValueError(f"Could not find {description}. Tried:\n{tried}")

def infer_project_root(config):
    """
    get OrthoSimulator.py project root
    """
    script_dir = os.path.dirname(os.path.abspath(__file__))
    return os.path.abspath(os.path.join(script_dir, "..", ".."))


def write_runtime_parameter_file(config, output_dir, project_root, threads):
    """
    write the runtime parameter file which LoadParameters reads
    """
    runtime_config = dict(config)
    runtime_config["threads"] = str(threads)

    runtime_config.setdefault(
        "iqtree_path",
        os.path.join(project_root, "bin", "iqtree3"),
    )
    runtime_config.setdefault(
        "sagephy_path",
        os.path.join(project_root, "bin", "sagephy-1.0.0.jar"),
    )
    runtime_config.setdefault("sigma_log_mean", "0")
    if "sigma_log_sd" not in runtime_config:
        if "max_sigma2" in runtime_config:
            runtime_config["sigma_log_sd"] = runtime_config["max_sigma2"]
        else:
            runtime_config["sigma_log_sd"] = "1"
    # ILS is off unless the user sets --ils-theta. When it is not set, no
    # ils_theta is written to the parameter file and the coalescent
    # rearrangement step (coalescence.py) is skipped entirely.

    runtime_parameter_file = os.path.join(
        output_dir,
        "orthogroup_simulation_parameters.txt",
    )
    with open(runtime_parameter_file, "w") as out:
        for key in sorted(runtime_config):
            out.write(f"{key}={runtime_config[key]}\n")

    return runtime_parameter_file

def configure_legacy_args(args, config):
    """
    translates the wrapper args into the things that this script needs
        args.o = absolute output directory
        args.n = number of orthogroups to simulate (int)
        args.s = path to the ultrametric species tree
        args.gap_profile = path to the empirical gap-position profile
        args.domain_species = species used for PFAM domain lookups
        args.p = runtime parameter file
        args.f = PFAM genome FASTA for domain_species
        args.d = PFAM domain-model CSV for domain_species
        args.PFAM = PFAM scan results .txt for domain_species
        args.project_root = OrthoSimulator.py project root
        args.pfam_root = root of the PFAM
    """
    project_root = infer_project_root(config)

    # core simulation inputs
    args.o = os.path.abspath(args.output)
    args.n = int(require_config(config, "Orthogroups"))
    args.s = require_config(config, "ultrametric_tree")
    args.gap_profile = require_config(config, "gap_file")
    domain_species = require_config(config, "species")
    args.domain_species = domain_species
    # Writes orthogroup_simulation_parameters.txt and returns its path.
    args.p = write_runtime_parameter_file(
        config=config,
        output_dir=args.o,
        project_root=project_root,
        threads=args.threads,
    )

    # locate PFAM data for domain_species
    pfam_root = config.get("pfam_dir")
    if pfam_root:
        pfam_root = resolve_path(pfam_root, project_root)
    else:
        pfam_root = os.path.join(project_root, "PFAM")

    # the PFAM data is shipped as one compressed archive per species;
    # unpack this species if its plain files are not there yet
    if pfam_archive.unpack(pfam_root, domain_species):
        print(f"  Unpacked PFAM data for {domain_species} into {pfam_root}")

    genome_candidates = []
    if config.get("pfam_genome_file"):
        genome_candidates.append(resolve_path(config["pfam_genome_file"], project_root))
    genome_candidates.extend([
        os.path.join(pfam_root, "pfam_genomes", f"{domain_species}.fa"),
        os.path.join(pfam_root, "pfam_genomes", f"{domain_species}.fasta"),
        os.path.join(pfam_root, "pfam_genomes", f"{domain_species}.faa"),
    ])
    args.f = first_existing_file(
        genome_candidates,
        f"PFAM genome FASTA for species '{domain_species}'",
    )

    domain_csv_candidates = []
    if config.get("pfam_domain_csv"):
        domain_csv_candidates.append(resolve_path(config["pfam_domain_csv"], project_root))
    domain_csv_candidates.extend([
        os.path.join(
            pfam_root,
            "pfam_results",
            "csv_files",
            "species_models_fixed",
            f"{domain_species}_pfam_models.csv",
        ),
        os.path.join(
            pfam_root,
            "pfam_results",
            "csv_files",
            "species_models_fixed",
            f"{domain_species}.csv",
        ),
    ])
    args.d = first_existing_file(
        domain_csv_candidates,
        f"PFAM domain-model CSV for species '{domain_species}'",
    )

    pfam_scan_candidates = []
    if config.get("pfam_scan_file"):
        pfam_scan_candidates.append(resolve_path(config["pfam_scan_file"], project_root))
    pfam_scan_candidates.extend([
        os.path.join(pfam_root, "pfam_results", f"{domain_species}_pfam.txt"),
        os.path.join(pfam_root, "pfam_results", f"{domain_species}.pfam.txt"),
        os.path.join(pfam_root, "pfam_results", f"{domain_species}.txt"),
        os.path.join(pfam_root, "pfam_results", "txt_files", f"{domain_species}_pfam.txt"),
        os.path.join(pfam_root, "pfamscan_results", f"{domain_species}_pfam.txt"),
    ])
    args.PFAM = first_existing_file(
        pfam_scan_candidates,
        f"PFAM scan results text file for species '{domain_species}'",
    )
    args.project_root = project_root
    args.pfam_root = pfam_root
    return args

def create_subfolders(output_dir):
    os.makedirs(output_dir, exist_ok=True)
    for folder in [
        "temporary_files",
        "proteome_files",
        "alignment_files",
        "tree_files",
    ]:
        os.makedirs(os.path.join(output_dir, folder), exist_ok=True)


def check_inputs_and_load_parameters(args):
    required_files = {
        "OrthoSimulator config": args.config,
        "species tree": args.s,
        "gap profile": args.gap_profile,
        "runtime parameter file": args.p,
        "PFAM genome FASTA": args.f,
        "PFAM scan results": args.PFAM,
        "PFAM domain-model CSV": args.d,
    }
    for description, path in required_files.items():
        if not os.path.isfile(path):
            raise ValueError(f"{description} not found: {path}")

    # Load runtime parameters and merge onto args.
    params = og.LoadParameters(args.p)
    for key, value in params.items():
        setattr(args, key, value)

    if not os.path.isfile(str(args.sagephy_path)):
        raise ValueError(f"SagePhy jar not found: {args.sagephy_path}")
    if not os.path.isfile(str(args.iqtree_path)):
        raise ValueError(f"IQ-TREE executable not found: {args.iqtree_path}")
    if not os.access(str(args.iqtree_path), os.X_OK):
        raise ValueError(f"IQ-TREE exists but is not executable: {args.iqtree_path}")

    # SaGePhy needs an ultrametric species tree (all tips the same distance
    # from the root) and otherwise stops with a Java stack trace
    check_species_tree_is_ultrametric(args.s)

    # SaGePhy rejects a birth-bias (gbc) of 0 or less and writes no trees
    gbc = getattr(args, "gbc", None)
    if gbc is not None and float(gbc) <= 0:
        raise ValueError(f"gbc must be > 0 (got {gbc}).")

    # Fail fast on a bad ils_theta rather than discovering it deep inside a
    # multiprocessing worker. Only checked when ILS was requested.
    ils_theta = getattr(args, "ils_theta", None)
    if ils_theta is not None and float(ils_theta) <= 0:
        raise ValueError(
            f"ils_theta must be > 0 (got {ils_theta}). A value of 0 is the "
            "extreme case of maximal ILS, not an off switch. To switch ILS "
            "off, leave out --ils-theta."
        )


def check_species_tree_is_ultrametric(tree_path, tolerance=1e-3):
    """
    Raise a clear error if the tips of the species tree are not all the same
    distance from the root. Rounding in the file is tolerated.
    """
    from ete4 import Tree
    tree = Tree(open(tree_path).read().strip(), parser=1)
    depths = []
    for leaf in tree.leaves():
        depth, node = 0.0, leaf
        while node.up is not None:
            depth += node.dist or 0.0
            node = node.up
        depths.append(depth)
    longest = max(depths)
    if longest <= 0:
        raise ValueError(f"The species tree has no branch lengths: {tree_path}")
    if (longest - min(depths)) / longest > tolerance:
        raise ValueError(
            f"The species tree is not ultrametric: tip depths range from {min(depths):.4g} to "
            f"{longest:.4g} in {tree_path}. OrthoSimulator needs a tree whose tips are all the "
            "same distance from the root."
        )


def print_detected_inputs(args):
    print("Detected OrthoSimulator-controlled inputs:")
    print(f"  Config:             {args.config}")
    print(f"  Project root:       {args.project_root}")
    print(f"  Output directory:   {args.o}")
    print(f"  Runtime parameters: {args.p}")
    print(f"  Orthogroups:        {args.n}")
    print(f"  Threads:            {args.threads}")
    print(f"  Species tree:       {args.s}")
    print(f"  Gap profile:        {args.gap_profile}")
    print(f"  PFAM root:          {args.pfam_root}")
    print(f"  PFAM species:       {args.domain_species}")
    print(f"  Genome FASTA:       {args.f}")
    print(f"  PFAM scan results:  {args.PFAM}")
    print(f"  Domain CSV:         {args.d}")
    print(f"  IQ-TREE:            {args.iqtree_path}")
    print(f"  SagePhy:            {args.sagephy_path}")


def run_orthogroup_simulation(
    complete_parameters,
    output_abolsute_path,
    threads,
    current_file_path,
    config_file_path,
):
    print("---------------------------------------------------")
    print("Welcome to the Orthogroup Simulation\n")

    ## greeting
    start_time = time.time()
    current_hour = time.localtime().tm_hour
    if 6 <= current_hour < 12:
        print("Good morning!")
    elif 12 <= current_hour < 18:
        print("Good afternoon!")
    else:
        print("Good evening!")

    ## threads
    threads = int(threads)
    if threads < 1:
        raise ValueError("threads must be >= 1")

    ## config and pathing
    config_file_path = os.path.abspath(config_file_path)
    output_abolsute_path = os.path.abspath(output_abolsute_path)

    if not os.path.isfile(config_file_path):
        raise ValueError(f"OrthoSimulator config not found: {config_file_path}")

    if not os.path.isdir(output_abolsute_path):
        raise ValueError(
            f"Output directory not found: {output_abolsute_path}\n"
            "OrthoSimulator.py should create this directory before running the simulation."
        )

    config = read_config(config_file_path)

    ##args
    args = SimpleNamespace(
        config=config_file_path,
        output=output_abolsute_path,
        threads=threads,
    )
    args = configure_legacy_args(args, config)

    create_subfolders(args.o)
    check_inputs_and_load_parameters(args)

    print_detected_inputs(args)

    ## multiprocesing
    mp.set_start_method("spawn", force=True)
    og.RunOrthogroupMultiProc(
        n=args.n,
        outdir=args.o,
        threads=threads,
        args=args,
    )

    # printing and saving output
    print("\n")

    skipped = og.SkippedOrthogroups(args)
    if skipped:
        print(f"WARNING: {len(skipped)} of {args.n} orthogroups were skipped because AliSim failed on them "
              f"(for example when a very deep tree deletes nearly every site). They are left out of every "
              f"output file. Skipped: {', '.join(skipped[:20])}{' ...' if len(skipped) > 20 else ''}\n")

    og.CopyAlignments(args)
    print("Written alignments\n")

    og.CopyTrees(args)
    print("Written trees\n")

    og.BuildProteomes(args)
    print("Written proteomes\n")

    og.ConcatOrthologs(args)
    print("Written orthologs\n")

    og.SaveOrthogroups(args)
    print("Written orthogroups\n")

    og.BuildOrthogroupStats(args)
    print("Written orthogroup stats\n")

    og.WriteRunLog(args)
    print("Written log file\n")

    elapsed = time.time() - start_time
    print(f"Output written to: {args.o}")
    print(f"Elapsed time: {elapsed:.2f} seconds\n")
    print("Have a smashing day!\n")
