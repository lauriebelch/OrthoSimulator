#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import csv
import re
from multiprocessing import Pool
#from ete3 import Tree
from ete4 import Tree

def main(orthofinder_folder, output_folder, n_threads):

    # --- DEFINE BASE PATHS ---
    input_base = os.path.abspath(orthofinder_folder)
    output_base = os.path.abspath(output_folder)

    # INPUT FILES
    orthogroups_file = os.path.join(
        input_base,
        "WorkingDirectory",
        "Numeric_Orthogroups.tsv",
    )

    proteome_dir = os.path.join(
        input_base,
        "WorkingDirectory",
    )

    # OUTPUT DIR
    output_dir = os.path.join(
        output_base,
        "TrainingResults",
        "Orthogroup_Sequences",
    )

    os.makedirs(output_dir, exist_ok=True)

    # --- STEP 1: build numeric orthogroups ---
    SpeciesDict, SequenceIDsDict = File_Dictionaries(input_base)

    # --- STEP 2: save matching numeric species tree ---
    Save_Numeric_Species_Tree(
        input_base=input_base,
        SpeciesDict=SpeciesDict,
    )

    # --- STEP 3: load numeric proteomes ---
    proteomes = {}

    for fname in os.listdir(proteome_dir):
        if re.fullmatch(r"Species\d+\.fa", fname):
            species = fname.split("Species")[1].split(".fa")[0]
            proteomes[species] = index_proteome_by_id(
                load_fasta(os.path.join(proteome_dir, fname))
            )

    # --- STEP 4: read numeric orthogroups ---
    with open(orthogroups_file) as f:
        reader = csv.DictReader(f, delimiter="\t")
        rows = list(reader)
        fieldnames = reader.fieldnames

    # --- STEP 5: write one sequence file per orthogroup ---
    with Pool(processes=int(n_threads)) as pool:
        pool.starmap(
            process_row,
            [
                (row, fieldnames, proteomes, output_dir)
                for row in rows
            ],
        )


# ============================================================
# Numeric species tree
# ============================================================

def Save_Numeric_Species_Tree(input_base, SpeciesDict):
    """
    Save a numeric species tree using the SAME species numbering used by
    Orthofinder's SpeciesIDs.txt / WorkingDirectory files.

    Input tree:
        orthofinder_folder/Species_Tree/SpeciesTree_rooted_node_labels.txt

    Outputs:
        orthofinder_folder/Species_Tree/SpeciesTree_rooted_node_labels_numeric.txt
        orthofinder_folder/Species_Tree/species_name_to_number.csv

    Example:
        species tree leaf: homo_sapiens
        SpeciesIDs code:   7
        numeric tree leaf: 7
    """

    species_tree_dir = os.path.join(
        input_base,
        "Species_Tree",
    )

    species_tree_src = os.path.join(
        species_tree_dir,
        "SpeciesTree_rooted_node_labels.txt",
    )

    numeric_species_tree_path = os.path.join(
        species_tree_dir,
        "SpeciesTree_rooted_node_labels_numeric.txt",
    )

    species_mapping_path = os.path.join(
        species_tree_dir,
        "species_name_to_number.csv",
    )

    if not os.path.isfile(species_tree_src):
        raise ValueError(
            f"Species tree not found: {species_tree_src}"
        )

    tree = Tree(species_tree_src, parser=1)

    for leaf in tree.leaves():
        leaf_name = leaf.name
        leaf_base = os.path.splitext(leaf_name)[0]

        if leaf_name in SpeciesDict:
            leaf.name = str(SpeciesDict[leaf_name])

        elif leaf_base in SpeciesDict:
            leaf.name = str(SpeciesDict[leaf_base])

        else:
            raise ValueError(
                f"Species '{leaf.name}' is in the species tree but was not "
                f"found in SpeciesIDs.txt"
            )

    tree.write(
        outfile=numeric_species_tree_path,
        parser=1,
    )

    with open(species_mapping_path, "w", newline="") as out:
        writer = csv.writer(out)
        writer.writerow(["species", "number"])

        for species in sorted(SpeciesDict):
            writer.writerow([
                species,
                SpeciesDict[species],
            ])

    print(f"Wrote numeric species tree to {numeric_species_tree_path}")
    print(f"Wrote species mapping to {species_mapping_path}")


# ============================================================
# Build numeric orthogroups
# ============================================================

def File_Dictionaries(base):

    OG_Path = os.path.join(
        base,
        "Orthogroups",
        "Orthogroups.tsv",
    )

    SeqIDs = os.path.join(
        base,
        "WorkingDirectory",
        "SequenceIDs.txt",
    )

    SpeciesIDs = os.path.join(
        base,
        "WorkingDirectory",
        "SpeciesIDs.txt",
    )

    Output = os.path.join(
        base,
        "WorkingDirectory",
        "Numeric_Orthogroups.tsv",
    )

    if os.path.exists(Output):
        os.remove(Output)

    SpeciesDict = {}
    Alt_SpeciesDict = {}

    with open(SpeciesIDs) as Species:
        for line in Species:
            if ": " not in line:
                continue

            key, _, value = line.partition(": ")

            raw_species = value.strip()
            species_base = os.path.splitext(raw_species)[0]
            species_code = key.strip()

            SpeciesDict[species_base] = species_code
            Alt_SpeciesDict[species_code] = species_base

    SequenceIDsDict = {
        code: {}
        for code in Alt_SpeciesDict
    }

    with open(SeqIDs) as SeqID:
        for line in SeqID:
            if ":" not in line:
                continue

            key, _, value = line.partition(":")

            sp_code = key.split("_")[0]
            coded_gene = key.strip()

            first_token = value.strip().split()[0]
            original_gene = first_token

            if sp_code in SequenceIDsDict:
                SequenceIDsDict[sp_code][original_gene] = coded_gene

    with open(OG_Path) as OG_file, open(Output, "w") as outfile:
        header = next(OG_file).rstrip("\n")
        colnames = header.split("\t")[1:]

        numeric_cols = [
            SpeciesDict[s]
            for s in colnames
        ]

        outfile.write(
            "Orthogroup\t" + "\t".join(numeric_cols) + "\n"
        )

        species_order = numeric_cols[:]

        for line in OG_file:
            if not line.startswith("OG"):
                outfile.write(line)
                continue

            parts = line.rstrip("\n").split("\t")

            og_name = parts[0]
            species_fields = parts[1:]
            new_fields = []

            for pos, field in enumerate(species_fields):
                if field == "":
                    new_fields.append("")
                    continue

                species_code = species_order[pos]

                genes = [
                    g for g in field.split(", ")
                    if g != ""
                ]

                replaced_genes = []

                for gene in genes:
                    try:
                        new_gene = SequenceIDsDict[species_code][gene]
                        replaced_genes.append(new_gene)

                    except KeyError:
                        raise KeyError(
                            f"Gene '{gene}' not found in SequenceIDs for "
                            f"species code '{species_code}'.\n"
                            f"Column species: {colnames[pos]}"
                        )

                new_fields.append(", ".join(replaced_genes))

            outfile.write(
                og_name + "\t" + "\t".join(new_fields) + "\n"
            )

    return SpeciesDict, SequenceIDsDict


# ============================================================
# FASTA loading
# ============================================================

def load_fasta(filepath):
    sequences = {}

    with open(filepath) as f:
        header = None
        seq = []

        for line in f:
            line = line.strip()

            if line.startswith(">"):
                if header:
                    sequences[header] = "".join(seq)

                header = line[1:]
                seq = []

            else:
                seq.append(line)

        if header:
            sequences[header] = "".join(seq)

    return sequences


def index_proteome_by_id(proteome_dict):
    """
    Build an exact-match index from a loaded proteome.

    Keys are the gene ID (the first whitespace-separated token of the
    FASTA header, matching how coded gene IDs are derived from
    SequenceIDs.txt elsewhere in this script). Using an exact key instead
    of substring matching avoids gene IDs like "0_5" incorrectly matching
    headers such as "10_50" or "0_51".

    Returns:
        {gene_id: (header, seq), ...}
    """
    indexed = {}

    for header, seq in proteome_dict.items():
        gene_id = header.split()[0]

        if gene_id in indexed:
            raise ValueError(
                f"Duplicate gene ID '{gene_id}' found in proteome "
                f"(headers '{indexed[gene_id][0]}' and '{header}')."
            )

        indexed[gene_id] = (header, seq)

    return indexed


def find_sequence(gene_id, proteome_index):
    return proteome_index.get(gene_id, (None, None))


# ============================================================
# Multiprocessing worker
# ============================================================

def process_row(row, fieldnames, proteomes, output_dir):
    og = row["Orthogroup"]

    output_path = os.path.join(
        output_dir,
        f"{og}.faa",
    )

    with open(output_path, "w") as out:
        for species in fieldnames[1:]:
            if not row[species]:
                continue

            genes = [
                g.strip()
                for g in row[species].split(",")
            ]

            for gene in genes:
                if not gene:
                    continue

                if species not in proteomes:
                    print(f"Warning: no proteome for {species}")
                    continue

                header, seq = find_sequence(
                    gene,
                    proteomes[species],
                )

                if header:
                    out.write(f">{header}\n{seq}\n")
                else:
                    print(f"Warning: {gene} not found in {species}")