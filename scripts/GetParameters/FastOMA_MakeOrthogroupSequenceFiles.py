#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import csv
import shutil
from collections import defaultdict

from Bio import SeqIO
#from ete3 import Tree
from ete4 import Tree


def main(
    fastoma_output,
    fastoma_input,
    output_base,
    n_threads,
    tool_tree=None,
):
    """
    Create numeric orthogroup sequence files and a matching numeric species tree.

    Important:
        The numeric species IDs used in gene names are the SAME numeric IDs
        used in the numeric species tree.

    Gene names become:
        1_originalGeneID
        2_originalGeneID
        etc.

    Numeric species tree leaves become:
        1, 2, 3, ...

    Parameters
    ----------
    fastoma_output : str
        Path to the FastOMA output/results folder.

    fastoma_input : str
        Path to the FastOMA input folder.

        Expected structure:
            fastoma_input/
            └── proteome/
                ├── Species1.fa
                ├── Species2.fa
                └── ...

    output_base : str
        OrthoSimulator output folder.

    n_threads : int
        Number of threads. Currently accepted for interface consistency.

    tool_tree : str or None
        Optional species tree path passed by OrthoSimulator.py via --tool-tree.

        If None, this script falls back to:
            fastoma_output/species_tree_checked.nwk
    """

    # ============================================================
    # Paths
    # ============================================================

    fastoma_out = os.path.abspath(fastoma_output)
    fastoma_in = os.path.abspath(fastoma_input)
    output_base = os.path.abspath(output_base)

    proteome_dir = os.path.join(
        fastoma_in,
        "proteome",
    )

    mapping_file = os.path.join(
        fastoma_out,
        "RootHOGs.tsv",
    )

    seq_dir = os.path.join(
        output_base,
        "TrainingResults",
        "Orthogroup_Sequences",
    )

    if tool_tree is None:
        species_tree_src = os.path.join(
            fastoma_out,
            "species_tree_checked.nwk",
        )
    else:
        species_tree_src = os.path.abspath(tool_tree)

    species_tree_dir = os.path.join(
        fastoma_out,
        "Species_Tree",
    )

    species_tree_dst = os.path.join(
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

    os.makedirs(seq_dir, exist_ok=True)
    os.makedirs(species_tree_dir, exist_ok=True)

    print(f"Mapping file: {mapping_file}")
    print(f"Proteomes: {proteome_dir}")
    print(f"Sequence output: {seq_dir}")
    print(f"Species tree: {species_tree_src}")

    # ============================================================
    # Validation
    # ============================================================

    if not os.path.isdir(fastoma_out):
        raise ValueError(f"FastOMA output directory not found: {fastoma_out}")

    if not os.path.isdir(fastoma_in):
        raise ValueError(f"FastOMA input directory not found: {fastoma_in}")

    if not os.path.isfile(mapping_file):
        raise ValueError(f"RootHOGs.tsv not found: {mapping_file}")

    if not os.path.isdir(proteome_dir):
        raise ValueError(f"Proteome directory not found: {proteome_dir}")

    if not os.path.isfile(species_tree_src):
        raise ValueError(f"Species tree not found: {species_tree_src}")

    # ============================================================
    # Load orthogroups
    # ============================================================

    orthogroups = defaultdict(list)

    with open(mapping_file) as f:
        next(f)

        for line in f:
            if not line.strip():
                continue

            og_id, gene_id = line.strip().split("\t")[:2]
            orthogroups[og_id].append(gene_id)

    print(f"Loaded {len(orthogroups)} orthogroups")

    # ============================================================
    # Index proteome sequences and create numeric gene IDs
    # ============================================================

    gene_numeric = {}

    species_to_id = {}
    species_counter = 0

    proteome_files = sorted(
        os.listdir(proteome_dir)
    )

    for proteome_file in proteome_files:
        if not proteome_file.endswith((".fa", ".faa", ".fasta")):
            continue

        species_name = os.path.splitext(proteome_file)[0]

        if species_name not in species_to_id:
            species_counter += 1
            species_to_id[species_name] = species_counter

        species_id = species_to_id[species_name]

        file_path = os.path.join(
            proteome_dir,
            proteome_file,
        )

        print(f"Indexing {species_name} as species {species_id}")

        for record in SeqIO.parse(file_path, "fasta"):
            original_id = record.id.split()[0]

            if original_id in gene_numeric:
                raise ValueError(
                    f"Duplicate gene ID '{original_id}' found while indexing "
                    f"proteome '{species_name}'. gene_numeric is keyed only "
                    f"by original ID (not species-qualified), so this gene "
                    f"would silently overwrite a previously indexed sequence "
                    f"and corrupt orthogroup sequence output. Check whether "
                    f"'{original_id}' is reused across multiple species "
                    f"proteomes."
                )

            rec_num = record[:]
            rec_num.id = f"{species_id}_{original_id}"
            rec_num.name = rec_num.id
            rec_num.description = ""

            gene_numeric[original_id] = rec_num

    print(f"Indexed {len(gene_numeric)} numeric sequences")
    print(f"Indexed {len(species_to_id)} species")

    if len(gene_numeric) == 0:
        raise ValueError(
            f"No protein sequences were indexed from: {proteome_dir}"
        )

    if len(species_to_id) == 0:
        raise ValueError(
            f"No species FASTA files were found in: {proteome_dir}"
        )

    # ============================================================
    # Write orthogroup sequence files
    # ============================================================

    written = 0
    missing_genes = 0

    for og_id, gene_ids in orthogroups.items():
        clean_id = og_id.replace("HOG:", "")

        if clean_id.startswith("OG"):
            og_name = clean_id
        else:
            og_name = f"OG{clean_id}"

        records_num = []

        for gene_id in gene_ids:
            if gene_id in gene_numeric:
                records_num.append(gene_numeric[gene_id])
            else:
                missing_genes += 1

        if len(records_num) < 2:
            continue

        out_fasta = os.path.join(
            seq_dir,
            f"{og_name}.fa",
        )

        SeqIO.write(
            records_num,
            out_fasta,
            "fasta",
        )

        written += 1

    print(f"Wrote {written} orthogroup sequence files")

    if missing_genes > 0:
        print(
            f"WARNING: {missing_genes} RootHOG gene IDs were not found "
            "in the indexed proteomes."
        )

    if written == 0:
        raise ValueError(
            "No orthogroup sequence files were written. "
            "Check that RootHOGs.tsv gene IDs match the proteome FASTA headers."
        )

    # ============================================================
    # Copy original species tree to compatible location
    # ============================================================

    shutil.copyfile(
        species_tree_src,
        species_tree_dst,
    )

    print(f"Copied species tree to {species_tree_dst}")

    # ============================================================
    # Save numeric species tree using SAME species_to_id mapping
    # ============================================================

    save_numeric_species_tree(
        species_tree_src=species_tree_src,
        numeric_species_tree_path=numeric_species_tree_path,
        species_mapping_path=species_mapping_path,
        species_to_id=species_to_id,
    )

    print("Sequence generation complete.")


def save_numeric_species_tree(
    species_tree_src,
    numeric_species_tree_path,
    species_mapping_path,
    species_to_id,
):
    """
    Save a copy of the species tree with leaf names replaced by the same
    numeric species IDs used in the orthogroup sequence FASTA files.

    Outputs:
        SpeciesTree_rooted_node_labels_numeric.txt
        species_name_to_number.csv
    """

    tree = Tree(
        species_tree_src,
        parser=1,
    )

    for leaf in tree.leaves():
        if leaf.name not in species_to_id:
            raise ValueError(
                f"Species '{leaf.name}' is in the species tree but was not "
                f"found in the proteome directory."
            )

        leaf.name = str(
            species_to_id[leaf.name]
        )

    tree.write(
        outfile=numeric_species_tree_path,
        parser=1,
    )

    with open(species_mapping_path, "w", newline="") as out:
        writer = csv.writer(out)

        writer.writerow(
            [
                "species",
                "number",
            ]
        )

        for species in sorted(species_to_id):
            writer.writerow(
                [
                    species,
                    species_to_id[species],
                ]
            )

    print(f"Wrote numeric species tree to {numeric_species_tree_path}")
    print(f"Wrote species mapping to {species_mapping_path}")