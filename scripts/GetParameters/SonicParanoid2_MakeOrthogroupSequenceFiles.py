#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import csv
import glob
from collections import defaultdict
from Bio import SeqIO
#from ete3 import Tree
from ete4 import Tree
import shutil


def main(sonic_output, sonic_input, output_base, n_threads, tool_tree):
    """
    Create numeric orthogroup sequence files and a matching numeric species tree
    from SonicParanoid2 output.

    tool_tree: path to the user-supplied species tree (required — SonicParanoid2
               does not produce a species tree natively).

    Expected input structure:

        sonic_output/
            runs/
                <run_name>/
                    ortholog_groups/
                        flat.ortholog_groups.tsv

        sonic_input/
            Species_A.fa
            Species_B.fa
            Species_C.fa
            ...

    SonicParanoid2 flat.ortholog_groups.tsv format:

        group_id    Species_A.fa    Species_B.fa    Species_C.fa
        34          gene1,gene2     *               gene3,gene4

    Output:

        output_base/TrainingResults/Orthogroup_Sequences/OG34.fa

    Gene names become:

        1_originalGeneID
        2_originalGeneID
        etc.

    Numeric species tree leaves become:

        1, 2, 3, ...

    Important:
        The numeric species IDs used in gene names are the SAME numeric IDs
        used in the numeric species tree.
    """

    # ============================================================
    # Paths
    # ============================================================

    sonic_out = os.path.abspath(sonic_output)
    sonic_in = os.path.abspath(sonic_input)
    output_base = os.path.abspath(output_base)

    proteome_dir = sonic_in

    groups_file = find_sonicparanoid_groups_file(sonic_out)

    species_tree_src = os.path.abspath(tool_tree)

    seq_dir = os.path.join(
        output_base,
        "TrainingResults",
        "Orthogroup_Sequences",
    )

    species_tree_dir = os.path.join(
        sonic_out,
        "Species_Tree",
    )

    species_tree_dst = os.path.join(
        species_tree_dir,
        "SpeciesTree_rooted_node_labels.txt",
    )

    numeric_species_tree_path = os.path.join(
    species_tree_dir,
    "SpeciesTree_rooted_node_labels_numeric.txt")
    
    species_mapping_path = os.path.join(
        species_tree_dir,
        "species_name_to_number.csv",
    )

    os.makedirs(seq_dir, exist_ok=True)
    os.makedirs(species_tree_dir, exist_ok=True)

    print(f"SonicParanoid2 output: {sonic_out}")
    print(f"SonicParanoid2 input proteomes: {proteome_dir}")
    print(f"Orthogroup file: {groups_file}")
    print(f"Species tree: {species_tree_src}")
    print(f"Sequence output: {seq_dir}")
    print(f"Species tree output: {species_tree_dst}")

    # ============================================================
    # Validation
    # ============================================================

    if not os.path.isdir(sonic_out):
        raise ValueError(f"SonicParanoid2 output directory not found: {sonic_out}")

    if not os.path.isdir(proteome_dir):
        raise ValueError(f"Proteome directory not found: {proteome_dir}")

    if not os.path.isfile(groups_file):
        raise ValueError(f"Orthogroup file not found: {groups_file}")

    if not os.path.isfile(species_tree_src):
        raise ValueError(f"Species tree not found: {species_tree_src}")

    # ============================================================
    # Load SonicParanoid2 orthogroups
    # ============================================================

    orthogroups, species_from_groups = load_sonicparanoid_orthogroups(groups_file)

    print(f"Loaded {len(orthogroups)} orthogroups")
    print(f"Found {len(species_from_groups)} species columns in orthogroup file")

    # ============================================================
    # Index proteomes and create numeric gene IDs
    # ============================================================

    gene_numeric, species_to_id = index_proteomes_numeric(
        proteome_dir=proteome_dir,
    )

    print(f"Indexed {len(gene_numeric)} numeric sequences")
    print(f"Indexed {len(species_to_id)} species")

    # ============================================================
    # Check species consistency
    # ============================================================

    check_species_consistency(
        species_from_groups=species_from_groups,
        species_to_id=species_to_id,
        species_tree_src=species_tree_src,
    )

    # ============================================================
    # Write orthogroup sequence files
    # ============================================================

    written = write_orthogroup_fastas(
        orthogroups=orthogroups,
        gene_numeric=gene_numeric,
        seq_dir=seq_dir,
    )

    print(f"Wrote {written} orthogroup sequence files")

    # ============================================================
    # Save numeric species tree using SAME species_to_id mapping
    # ============================================================
    shutil.copy2(species_tree_src, species_tree_dst)

    save_numeric_species_tree(
        species_tree_src=species_tree_src,
        numeric_species_tree_path=numeric_species_tree_path,
        species_mapping_path=species_mapping_path,
        species_to_id=species_to_id,
    )

    print("SonicParanoid2 sequence generation complete.")


def find_sonicparanoid_groups_file(sonic_output):
    """
    Find SonicParanoid2 flat.ortholog_groups.tsv.

    Expected pattern:

        sonic_output/runs/*/ortholog_groups/flat.ortholog_groups.tsv

    If multiple files are found, this function raises an error so that you do
    not accidentally use the wrong SonicParanoid2 run.
    """

    pattern = os.path.join(
        sonic_output,
        "runs",
        "*",
        "ortholog_groups",
        "flat.ortholog_groups.tsv",
    )

    matches = sorted(glob.glob(pattern))

    if len(matches) == 0:
        raise ValueError(
            f"No flat.ortholog_groups.tsv file found with pattern:\n{pattern}"
        )

    if len(matches) > 1:
        raise ValueError(
            "Multiple flat.ortholog_groups.tsv files found:\n"
            + "\n".join(matches)
            + "\nPlease keep only one run in sonic_output/runs/, "
              "or modify the script to select a specific run."
        )

    return matches[0]


def load_sonicparanoid_orthogroups(groups_file):
    """
    Read SonicParanoid2 flat.ortholog_groups.tsv.

    Expected format:

        group_id    Species_A.fa    Species_B.fa    Species_C.fa
        34          gene1,gene2     *               gene3,gene4

    Missing entries are represented by '*'.

    Returns:
        orthogroups:
            dict:
                group_id -> [(species_name, gene_id), ...]

        species_from_groups:
            set of species names found in the orthogroup file columns

    Notes:
        Species names are taken from column names by removing the file extension.

        Example:
            Homo_sapiens.fa -> Homo_sapiens
    """

    orthogroups = defaultdict(list)
    species_from_groups = set()

    with open(groups_file) as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader)

        if not header:
            raise ValueError(f"Empty header in {groups_file}")

        if header[0] != "group_id":
            raise ValueError(
                f"Expected first column to be 'group_id', got '{header[0]}'"
            )

        species_columns = header[1:]

        species_names = [
            os.path.splitext(col)[0]
            for col in species_columns
        ]

        species_from_groups.update(species_names)

        for row in reader:
            if not row:
                continue

            group_id = row[0].strip()

            if not group_id:
                continue

            # If a row is shorter than the header, pad it with empty cells.
            # This avoids zip silently dropping missing trailing species.
            if len(row) < len(header):
                row = row + [""] * (len(header) - len(row))

            for species_name, cell in zip(species_names, row[1:]):
                cell = cell.strip()

                if not cell or cell == "*":
                    continue

                gene_ids = [
                    g.strip()
                    for g in cell.split(",")
                    if g.strip() and g.strip() != "*"
                ]

                for gene_id in gene_ids:
                    orthogroups[group_id].append(
                        (
                            species_name,
                            gene_id,
                        )
                    )

    return orthogroups, species_from_groups


def index_proteomes_numeric(proteome_dir):
    """
    Index input proteomes and create numeric gene IDs.

    Returns:
        gene_numeric:
            dict:
                (species_name, original_gene_id) -> numeric SeqRecord

        species_to_id:
            dict:
                species_name -> numeric species ID

    Gene IDs are species-aware, so the same gene ID can safely occur in
    different species without collisions.
    """

    gene_numeric = {}

    species_to_id = {}
    species_counter = 0

    proteome_files = sorted(os.listdir(proteome_dir))

    for proteome_file in proteome_files:
        if not proteome_file.endswith((".fa", ".fasta", ".faa")):
            continue

        species_name = os.path.splitext(proteome_file)[0]

        if species_name in species_to_id:
            raise ValueError(
                f"Duplicate species name derived from proteome files: {species_name}"
            )

        species_counter += 1
        species_to_id[species_name] = species_counter

        species_id = species_to_id[species_name]
        file_path = os.path.join(proteome_dir, proteome_file)

        print(f"Indexing {species_name} as species {species_id}")

        for record in SeqIO.parse(file_path, "fasta"):
            original_id = record.id.split()[0]

            key = (
                species_name,
                original_id,
            )

            if key in gene_numeric:
                raise ValueError(
                    f"Duplicate gene ID within species {species_name}: {original_id}"
                )

            rec_num = record[:]
            rec_num.id = f"{species_id}_{original_id}"
            rec_num.name = rec_num.id
            rec_num.description = ""

            gene_numeric[key] = rec_num

    if not species_to_id:
        raise ValueError(f"No proteome FASTA files found in {proteome_dir}")

    return gene_numeric, species_to_id


def check_species_consistency(
    species_from_groups,
    species_to_id,
    species_tree_src,
):
    """
    Check that species names are consistent among:

        1. SonicParanoid2 orthogroup file columns
        2. input proteome filenames
        3. species tree leaves

    This assumes all names are filename stems without .fa.
    """

    species_from_proteomes = set(species_to_id.keys())

    missing_from_proteomes = sorted(
        species_from_groups - species_from_proteomes
    )

    extra_in_proteomes = sorted(
        species_from_proteomes - species_from_groups
    )

    if missing_from_proteomes:
        raise ValueError(
            "These species occur in the SonicParanoid2 orthogroup file "
            "but no matching proteome file was found:\n"
            + "\n".join(missing_from_proteomes)
        )

    if extra_in_proteomes:
        print(
            "Warning: these proteome species do not occur in the "
            "SonicParanoid2 orthogroup file:\n"
            + "\n".join(extra_in_proteomes)
        )

    tree = Tree(species_tree_src, parser=1)

    species_from_tree = set(
        leaf.name
        for leaf in tree.leaves()
    )

    missing_from_tree = sorted(
        species_from_proteomes - species_from_tree
    )

    extra_in_tree = sorted(
        species_from_tree - species_from_proteomes
    )

    if missing_from_tree:
        raise ValueError(
            "These proteome species are missing from the species tree:\n"
            + "\n".join(missing_from_tree)
        )

    if extra_in_tree:
        raise ValueError(
            "These species are present in the species tree but not in "
            "the proteome directory:\n"
            + "\n".join(extra_in_tree)
        )


def write_orthogroup_fastas(
    orthogroups,
    gene_numeric,
    seq_dir,
):
    """
    Write one numeric FASTA file per SonicParanoid2 orthogroup.

    Orthogroup members are species-aware:

        (species_name, gene_id)

    Output filenames:

        OG34.fa
        OG35.fa
        ...
    """

    written = 0
    total_missing = 0
    groups_with_missing = 0

    for group_id, members in orthogroups.items():
        clean_id = str(group_id).strip()

        # SonicParanoid2 group IDs are usually numeric, so this gives OG34.
        # If a future group ID already starts with OG, avoid OGO...
        if clean_id.startswith("OG"):
            og_name = clean_id
        else:
            og_name = f"OG{clean_id}"

        records_num = []
        missing = []

        for species_name, gene_id in members:
            key = (
                species_name,
                gene_id,
            )

            if key in gene_numeric:
                records_num.append(gene_numeric[key])
            else:
                missing.append(key)

        if missing:
            groups_with_missing += 1
            total_missing += len(missing)

            print(
                f"Warning: {og_name} has {len(missing)} genes missing from proteomes"
            )

        if len(records_num) < 2:
            continue

        out_fasta = os.path.join(seq_dir, f"{og_name}.fa")

        SeqIO.write(
            records_num,
            out_fasta,
            "fasta",
        )

        written += 1

    if total_missing > 0:
        print(
            f"Warning: {total_missing} genes were missing across "
            f"{groups_with_missing} orthogroups"
        )

    return written


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
        SpeciesTree_rooted_node_labels.txt
        species_name_to_number.csv
    """

    tree = Tree(species_tree_src, parser=1)

    for leaf in tree.leaves():
        if leaf.name not in species_to_id:
            raise ValueError(
                f"Species '{leaf.name}' is in the species tree but was not "
                f"found in the proteome directory."
            )

        leaf.name = str(species_to_id[leaf.name])

    tree.write(
        outfile=numeric_species_tree_path,
        parser=1,
    )

    with open(species_mapping_path, "w", newline="") as out:
        writer = csv.writer(out)
        writer.writerow(["species", "number"])

        for species in sorted(species_to_id):
            writer.writerow([
                species,
                species_to_id[species],
            ])

    print(f"Wrote numeric species tree to {numeric_species_tree_path}")
    print(f"Wrote species mapping to {species_mapping_path}")