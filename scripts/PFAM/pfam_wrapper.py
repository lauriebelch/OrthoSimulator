#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Internal PFAM workflow wrapper for OrthoSimulator.

This script is called by OrthoSimulator.py, not normally by users directly.

Important layout:

    scripts/PFAM/
        pfam_wrapper.py
        PfamModelExtractor.py
        PfamCSVWriter.py
        download_pfamDB.sh
        run_pfamscan.sh

    PFAM/
        pfam_db/
        pfam_genomes/
        pfamscan_results/
        pfam_results/

So:
    scripts/PFAM = code only
    PFAM         = generated data/results only

Important behaviour:
    Only the FASTA files supplied in the current --Proteome call are processed.
    Existing files in PFAM/pfam_genomes or PFAM/pfamscan_results are not swept up
    into the current run.

Workflow:
    1. Stage current input FASTA file(s) into PFAM/pfam_genomes.
    2. Ensure PFAM/pfam_db exists and is hmmpress-ready.
    3. Run pfam_scan.pl only for the FASTA file(s) staged in this run.
    4. Copy this run's pfam_scan outputs into the OrthoSimulator run folder.
    5. Run PfamModelExtractor.py only on this run's scan outputs.
    6. Run PfamCSVWriter.py only on this run's species.
"""

import argparse
import importlib.util
import shutil
import subprocess
import sys
from pathlib import Path


FASTA_EXTENSIONS = {".fa", ".fasta", ".faa", ".fas", ".fna"}

PFAM_DB_REQUIRED_FILES = [
    "Pfam-A.hmm",
    "Pfam-A.hmm.dat",
    "Pfam-A.hmm.h3f",
    "Pfam-A.hmm.h3i",
    "Pfam-A.hmm.h3m",
    "Pfam-A.hmm.h3p",
]


def parse_args():
    parser = argparse.ArgumentParser(
        description="Internal OrthoSimulator PFAM workflow wrapper."
    )

    parser.add_argument(
        "--proteome",
        required=True,
        help="Input proteome FASTA file, or a directory containing FASTA files.",
    )

    parser.add_argument(
        "--threads",
        default="1",
        help="Number of CPU threads to pass to PFAM scan and IQ-TREE.",
    )

    parser.add_argument(
        "--pfam-dir",
        required=True,
        help=(
            "Top-level PFAM results/data directory. "
            "This should normally be OrthoSimulator/PFAM, not scripts/PFAM."
        ),
    )

    parser.add_argument(
        "--run-output",
        required=True,
        help=(
            "Outer OrthoSimulator output folder. Used to isolate the current run's "
            "PFAM scan inputs so old scans are not processed accidentally."
        ),
    )

    return parser.parse_args()


def safe_int_threads(value):
    try:
        threads = int(value)
    except ValueError:
        print(f"ERROR: --threads must be an integer, got: {value}")
        sys.exit(1)

    if threads < 1:
        print(f"ERROR: --threads must be >= 1, got: {value}")
        sys.exit(1)

    return threads


def run_command(cmd, description):
    print(f"\n{description}")
    print(" ".join(str(x) for x in cmd))

    try:
        subprocess.run(cmd, check=True)

    except FileNotFoundError as err:
        print(f"ERROR: command not found while running: {description}")
        print(err)
        sys.exit(1)

    except subprocess.CalledProcessError as err:
        print(f"ERROR: command failed while running: {description}")
        sys.exit(err.returncode)


def make_dirs(paths):
    for path in paths:
        path.mkdir(parents=True, exist_ok=True)


def is_fasta_file(path):
    return path.is_file() and path.suffix.lower() in FASTA_EXTENSIONS


def find_input_fastas(proteome_path):
    proteome_path = Path(proteome_path).expanduser().resolve()

    if not proteome_path.exists():
        print(f"ERROR: Proteome path does not exist: {proteome_path}")
        sys.exit(1)

    if proteome_path.is_file():
        if not is_fasta_file(proteome_path):
            print("ERROR: Proteome input is a file but does not look like a FASTA file:")
            print(proteome_path)
            print("Supported extensions: " + ", ".join(sorted(FASTA_EXTENSIONS)))
            sys.exit(1)

        return [proteome_path]

    fastas = sorted(
        path for path in proteome_path.iterdir()
        if is_fasta_file(path)
    )

    if not fastas:
        print("ERROR: no FASTA files found in proteome directory:")
        print(proteome_path)
        print("Supported extensions: " + ", ".join(sorted(FASTA_EXTENSIONS)))
        sys.exit(1)

    return fastas


def stage_fastas(input_fastas, pfam_genomes_dir):
    """
    Copy only the current input FASTA files into PFAM/pfam_genomes.

    Existing files with the same name are overwritten so repeated PFAM runs use
    the user's latest input. Other previously staged PFAM files are left in
    place but are not processed by this wrapper run.
    """
    staged = []

    for fasta in input_fastas:
        destination = pfam_genomes_dir / fasta.name

        if fasta.resolve() != destination.resolve():
            shutil.copy2(fasta, destination)

        staged.append(destination)
        print(f"Staged proteome for current run: {destination}")

    return staged


def species_name_from_scan_file(path):
    name = path.name
    if name.endswith("_pfam.txt"):
        return name.replace("_pfam.txt", "")
    return path.stem


def pfam_db_is_ready(pfam_db_dir):
    return all(
        (pfam_db_dir / filename).is_file()
        for filename in PFAM_DB_REQUIRED_FILES
    )


def ensure_pfam_db(script_dir, pfam_db_dir):
    """
    Ensure the PFAM database exists in PFAM/pfam_db.

    The download script lives in scripts/PFAM, so this function receives
    script_dir, not pfam_dir.
    """
    if pfam_db_is_ready(pfam_db_dir):
        print(f"PFAM database already prepared: {pfam_db_dir}")
        return

    print(f"PFAM database missing or incomplete: {pfam_db_dir}")

    script = script_dir / "download_pfamDB.sh"

    if not script.is_file():
        print("ERROR: PFAM database download script not found:")
        print(script)
        sys.exit(1)

    run_command(
        ["bash", str(script), str(pfam_db_dir)],
        "Downloading and preparing PFAM database",
    )

    if not pfam_db_is_ready(pfam_db_dir):
        print("ERROR: PFAM database setup completed, but required files are still missing.")
        print("Expected files:")

        for filename in PFAM_DB_REQUIRED_FILES:
            print(f"  {pfam_db_dir / filename}")

        sys.exit(1)


def prepare_current_run_scan_dir(run_output):
    output_dir = Path(run_output).expanduser().resolve()

    if not output_dir.is_dir():
        print("ERROR: OrthoSimulator output folder does not exist:")
        print(output_dir)
        sys.exit(1)

    current_run_pfamscan_dir = output_dir / "pfamscan_results_current_run"

    # This folder should describe exactly this invocation, so clear any leftovers
    # from an interrupted/repeated run using the same --output folder.
    if current_run_pfamscan_dir.exists():
        shutil.rmtree(current_run_pfamscan_dir)

    current_run_pfamscan_dir.mkdir(parents=True, exist_ok=True)
    return current_run_pfamscan_dir


def run_pfamscan_for_fastas(
    staged_fastas,
    pfam_db_dir,
    pfamscan_results_dir,
    current_run_pfamscan_dir,
    threads,
):
    """
    Run pfam_scan.pl only for FASTA files in staged_fastas.

    Each scan is written to the canonical PFAM/pfamscan_results folder, then
    copied into the current OrthoSimulator run folder. Downstream extraction uses the
    current-run folder, not the canonical folder, so old scan files are ignored.
    """
    current_scan_files = []

    for fasta in staged_fastas:
        outfile = pfamscan_results_dir / f"{fasta.stem}_pfam.txt"
        current_outfile = current_run_pfamscan_dir / outfile.name

        cmd = [
            "pfam_scan.pl",
            "-fasta", str(fasta),
            "-dir", str(pfam_db_dir),
            "-cpu", str(threads),
            "-outfile", str(outfile),
        ]

        run_command(cmd, f"Running pfam_scan.pl for current input {fasta.name}")

        if not outfile.is_file():
            print("ERROR: pfam_scan.pl finished but expected output was not created:")
            print(outfile)
            sys.exit(1)

        shutil.copy2(outfile, current_outfile)
        current_scan_files.append(current_outfile)
        print(f"Current-run PFAM scan copy: {current_outfile}")

    return current_scan_files


def run_model_extractor(script_dir, current_run_pfamscan_dir, pfam_results_dir, threads):
    """
    Run PfamModelExtractor directly (no subprocess).

    It receives only the current-run pfam_scan output folder, not the canonical
    PFAM/pfamscan_results folder that may contain old files.
    """
    script = script_dir / "PfamModelExtractor.py"

    if not script.is_file():
        print("ERROR: PFAM model extractor not found:")
        print(script)
        sys.exit(1)

    print("\nExtracting PFAM substitution models for current input only")

    spec = importlib.util.spec_from_file_location("PfamModelExtractor", str(script))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)

    extractor = mod.PfamModelExtractor(
        str(current_run_pfamscan_dir),
        str(pfam_results_dir),
        threads,
    )
    extractor.run()


def run_csv_writer(script_dir, current_run_pfamscan_dir, pfam_results_dir):
    """
    Run PfamCSVWriter directly (no subprocess).

    The CSV writer receives only the current-run scan folder, so it only writes
    gene-level model CSVs for species/genomes supplied in this run.
    """
    script = script_dir / "PfamCSVWriter.py"

    if not script.is_file():
        print("ERROR: PFAM CSV writer not found:")
        print(script)
        sys.exit(1)

    species_models_dir = pfam_results_dir / "csv_files" / "species_models"
    fixed_out_dir = pfam_results_dir / "csv_files" / "species_models_fixed"

    current_species = [
        species_name_from_scan_file(path)
        for path in sorted(current_run_pfamscan_dir.glob("*.txt"))
    ]

    print("\nWriting gene-level PFAM model CSVs for current input only")

    spec = importlib.util.spec_from_file_location("PfamCSVWriter", str(script))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)

    mod.run(
        pfam_scan_folder=str(current_run_pfamscan_dir),
        species_models_dir=str(species_models_dir),
        fixed_out_dir=str(fixed_out_dir),
        species=",".join(current_species),
    )


def write_run_record(
    run_output,
    pfam_dir,
    staged_fastas,
    current_run_pfamscan_dir,
    current_scan_files,
):
    output_dir = Path(run_output).expanduser().resolve()

    if not output_dir.is_dir():
        return

    record = output_dir / "pfam_run_paths.txt"

    with open(record, "w") as handle:
        handle.write(f"PFAM_DIR={pfam_dir}\n")
        handle.write(f"PFAM_GENOMES={pfam_dir / 'pfam_genomes'}\n")
        handle.write(f"PFAMSCAN_RESULTS={pfam_dir / 'pfamscan_results'}\n")
        handle.write(f"PFAMSCAN_RESULTS_CURRENT_RUN={current_run_pfamscan_dir}\n")
        handle.write(f"PFAM_RESULTS={pfam_dir / 'pfam_results'}\n")
        handle.write(f"PFAM_FINAL_CSV_DIR={pfam_dir / 'pfam_results' / 'csv_files' / 'species_models_fixed'}\n")
        handle.write("STAGED_FASTAS=\n")

        for fasta in staged_fastas:
            handle.write(f"  {fasta}\n")

        handle.write("CURRENT_RUN_PFAMSCAN_FILES=\n")

        for scan_file in current_scan_files:
            handle.write(f"  {scan_file}\n")

    print(f"Wrote PFAM run path record: {record}")


def run(proteome, threads, pfam_dir, run_output):
    """
    Direct callable entry point (used by OrthoSimulator.py without subprocess).

    Parameters
    ----------
    proteome : str
        Input proteome FASTA file or directory of FASTA files.
    threads : int or str
        Number of CPU threads.
    pfam_dir : str
        Top-level PFAM results/data directory (OrthoSimulator/PFAM).
    run_output : str
        OrthoSimulator run output folder.
    """
    threads = safe_int_threads(threads)

    # Code lives here:
    #   OrthoSimulator/scripts/PFAM
    script_dir = Path(__file__).resolve().parent

    # Results/data live here:
    #   OrthoSimulator/PFAM
    pfam_dir = Path(pfam_dir).expanduser().resolve()

    pfam_db_dir = pfam_dir / "pfam_db"
    pfam_genomes_dir = pfam_dir / "pfam_genomes"
    pfamscan_results_dir = pfam_dir / "pfamscan_results"
    pfam_results_dir = pfam_dir / "pfam_results"

    make_dirs([
        pfam_dir,
        pfam_db_dir,
        pfam_genomes_dir,
        pfamscan_results_dir,
        pfam_results_dir,
    ])

    current_run_pfamscan_dir = prepare_current_run_scan_dir(run_output)

    print("Starting OrthoSimulator PFAM workflow")
    print(f"PFAM scripts:                 {script_dir}")
    print(f"PFAM directory:               {pfam_dir}")
    print(f"PFAM database:                {pfam_db_dir}")
    print(f"PFAM genomes:                 {pfam_genomes_dir}")
    print(f"PFAM scan results:            {pfamscan_results_dir}")
    print(f"Current-run PFAM scan input:  {current_run_pfamscan_dir}")
    print(f"PFAM model results:           {pfam_results_dir}")
    print(f"Threads:                      {threads}")

    input_fastas = find_input_fastas(proteome)
    staged_fastas = stage_fastas(input_fastas, pfam_genomes_dir)

    ensure_pfam_db(script_dir, pfam_db_dir)

    current_scan_files = run_pfamscan_for_fastas(
        staged_fastas,
        pfam_db_dir,
        pfamscan_results_dir,
        current_run_pfamscan_dir,
        threads,
    )

    run_model_extractor(
        script_dir,
        current_run_pfamscan_dir,
        pfam_results_dir,
        threads,
    )

    run_csv_writer(
        script_dir,
        current_run_pfamscan_dir,
        pfam_results_dir,
    )

    write_run_record(
        run_output,
        pfam_dir,
        staged_fastas,
        current_run_pfamscan_dir,
        current_scan_files,
    )

    print("\nPFAM workflow complete.")
    print("Final downstream-ready CSV folder:")
    print(pfam_results_dir / "csv_files" / "species_models_fixed")
    print("Current-run PFAM scan folder:")
    print(current_run_pfamscan_dir)


def main():
    args = parse_args()
    run(args.proteome, args.threads, args.pfam_dir, args.run_output)


if __name__ == "__main__":
    main()
