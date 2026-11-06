#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import csv
import os
from pathlib import Path


def parse_args():
    script_dir = Path(__file__).resolve().parent
    # scripts/PFAM contains code only. Defaults point to the top-level PFAM
    # data/results folder when this script is run manually from the normal
    # OrthoSimulator layout. scripts/PFAM/pfam_wrapper.py still passes these explicitly.
    project_root = script_dir.parents[1]
    default_pfam_dir = project_root / "PFAM"
    default_pfam_scan = default_pfam_dir / "pfamscan_results"
    default_species_models = default_pfam_dir / "pfam_results" / "csv_files" / "species_models"
    default_fixed = default_pfam_dir / "pfam_results" / "csv_files" / "species_models_fixed"

    parser = argparse.ArgumentParser(
        description="Convert PFAM-level model CSVs to gene/file-level model CSVs."
    )
    parser.add_argument(
        "--pfam-scan-folder",
        default=str(default_pfam_scan),
        help="Folder containing pfam_scan.pl .txt output files.",
    )
    parser.add_argument(
        "--species-models-dir",
        default=str(default_species_models),
        help="Folder containing *_pfam_models.csv files from PfamModelExtractor.py.",
    )
    parser.add_argument(
        "--fixed-out-dir",
        default=str(default_fixed),
        help="Output folder for gene/file-level PFAM model CSVs.",
    )
    parser.add_argument(
        "--species",
        default=None,
        help=(
            "Optional comma-separated species/genome names to write. "
            "When omitted, species are inferred from --pfam-scan-folder."
        ),
    )
    return parser.parse_args()


def species_from_model_csv(filename):
    """
    Example:
        Mnemiopsis_leidyi_pfam_models.csv -> Mnemiopsis_leidyi
    """
    return filename.replace("_pfam_models.csv", "")


def species_from_pfam_scan_file(path):
    """
    Example:
        Mnemiopsis_leidyi_pfam.txt -> Mnemiopsis_leidyi
    """
    name = path.name
    if name.endswith("_pfam.txt"):
        return name.replace("_pfam.txt", "")
    return path.stem


def parse_species_list(species_arg, pfamscan_dir):
    if species_arg is not None:
        species = [item.strip() for item in species_arg.split(",") if item.strip()]
        return sorted(set(species))

    return sorted(
        species_from_pfam_scan_file(path)
        for path in pfamscan_dir.glob("*.txt")
    )


def pfamscan_file_for_species(species, pfamscan_dir):
    """
    Example:
        Mnemiopsis_leidyi -> Mnemiopsis_leidyi_pfam.txt
    """
    return pfamscan_dir / f"{species}_pfam.txt"


def build_pfam_to_genes_map(pfamscan_file):
    """
    Read pfam_scan.pl output.

    Column 1 = gene name
    Column 6 = PFAM accession, e.g. PF16021.10
    Column 8 = type, which should be Domain for the rows used here.
    """
    pfam_to_genes = {}

    with open(pfamscan_file) as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            parts = line.split()
            if len(parts) < 8:
                continue

            gene_name = parts[0]
            pfam_raw = parts[5]
            hit_type = parts[7]

            if hit_type != "Domain":
                continue

            pfam_id = pfam_raw.split(".")[0]
            pfam_to_genes.setdefault(pfam_id, [])
            if gene_name not in pfam_to_genes[pfam_id]:
                pfam_to_genes[pfam_id].append(gene_name)

    return pfam_to_genes


def fix_species_model_csv(model_csv, pfamscan_file, out_csv):
    """
    Replace first column of species model CSV.

    Input first column:
        pfam_id

    Output first column:
        file

    For each PFAM model row, write one row per gene containing that PFAM.
    """
    pfam_to_genes = build_pfam_to_genes_map(pfamscan_file)

    with open(model_csv) as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)
        old_fields = reader.fieldnames

    if old_fields is None:
        raise ValueError(f"No header found in {model_csv}")

    if old_fields[0] != "pfam_id":
        raise ValueError(
            f"Expected first column to be 'pfam_id' in {model_csv}, "
            f"but found '{old_fields[0]}'"
        )

    new_fields = ["file"] + old_fields[1:]
    fixed_rows = []

    for row in rows:
        pfam_id = row["pfam_id"]
        genes = pfam_to_genes.get(pfam_id, [])
        if not genes:
            continue

        for gene in genes:
            new_row = {"file": gene}
            for field in old_fields[1:]:
                new_row[field] = row[field]
            fixed_rows.append(new_row)

    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with open(out_csv, "w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=new_fields)
        writer.writeheader()
        writer.writerows(fixed_rows)

    return len(fixed_rows)


def run(pfam_scan_folder, species_models_dir, fixed_out_dir, species=None):
    """
    Direct callable entry point (used by pfam_wrapper.py without subprocess).

    Parameters
    ----------
    pfam_scan_folder : str or Path
        Folder containing pfam_scan.pl .txt output files.
    species_models_dir : str or Path
        Folder containing *_pfam_models.csv files from PfamModelExtractor.py.
    fixed_out_dir : str or Path
        Output folder for gene/file-level PFAM model CSVs.
    species : str or None
        Optional comma-separated species names. Auto-detected from scan folder
        when None.
    """
    pfamscan_dir       = Path(pfam_scan_folder).expanduser().resolve()
    species_models_dir = Path(species_models_dir).expanduser().resolve()
    fixed_out_dir      = Path(fixed_out_dir).expanduser().resolve()

    fixed_out_dir.mkdir(parents=True, exist_ok=True)

    current_species = parse_species_list(species, pfamscan_dir)

    print(f"Species model folder: {species_models_dir}")
    print(f"PFAM scan folder:     {pfamscan_dir}")
    print(f"Fixed output folder:  {fixed_out_dir}")
    print(f"Current-run species:  {len(current_species)}")

    if not current_species:
        print("WARNING: no current-run species found to convert.")
        print("Done.")
        return

    for sp in current_species:
        model_csv     = species_models_dir / f"{sp}_pfam_models.csv"
        pfamscan_file = pfamscan_file_for_species(sp, pfamscan_dir)
        out_csv       = fixed_out_dir / f"{sp}_pfam_models.csv"

        if not pfamscan_file.is_file():
            print(f"WARNING: no current-run PFAM scan file for {sp}: {pfamscan_file}")
            continue

        if not model_csv.is_file():
            print(f"WARNING: no species model CSV for {sp}: {model_csv}")
            continue

        n_rows = fix_species_model_csv(
            model_csv=model_csv,
            pfamscan_file=pfamscan_file,
            out_csv=out_csv,
        )
        print(f"{sp}: wrote {n_rows} rows")

    print("Done.")


def main():
    args = parse_args()
    run(
        pfam_scan_folder=args.pfam_scan_folder,
        species_models_dir=args.species_models_dir,
        fixed_out_dir=args.fixed_out_dir,
        species=args.species,
    )


if __name__ == "__main__":
    main()
