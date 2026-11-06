#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import os
import time
import sys
import shutil

# Ensure the script's directory is in PYTHONPATH before any local imports.
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
if SCRIPT_DIR not in sys.path:
    sys.path.insert(0, SCRIPT_DIR)

import MakeOrthogroupSequenceFiles
import CallDuplications
import EstimateParameters


def run(tool_output, output_folder_path, n_threads, tools_proteomes=None, tool_tree=None):
    ortho_folder_path  = os.path.abspath(tool_output)
    output_folder_path = os.path.abspath(output_folder_path)

    os.makedirs(output_folder_path, exist_ok=True)

    # --- Welcome ---
    print("---------------------------------------------------")
    print("Welcome to OrthoSimulator\n")
    print("Running GetParameters for OrthoFinder")
    print(f"OrthoFinder output: {ortho_folder_path}")
    print(f"OrthoSimulator output:  {output_folder_path}")
    print(f"Threads:            {n_threads}")
    print("")

    start_time = time.time()
    current_hour = time.localtime().tm_hour
    if 6 <= current_hour < 12:
        print("Good morning!")
    elif 12 <= current_hour < 18:
        print("Good afternoon!")
    else:
        print("Good evening!")
    print("---------------------------------------------------")

    # --- Pipeline ---
    print("Making Orthogroup sequence files ...")
    MakeOrthogroupSequenceFiles.main(
        ortho_folder_path,
        output_folder_path,
        n_threads,
    )
    print("Orthogroup sequence files created.")
    print("---------------------------------------------------")

    print("Aligning and Treeing to Estimate parameters.")
    done_file = os.path.join(output_folder_path, "TrainingResults", "Orthogroup_analysis", "DONE")
    if os.path.exists(done_file):
        print("Parameter estimation already completed — skipping.")
    else:
        EstimateParameters.main(
            ortho_folder_path,
            output_folder_path,
            n_threads,
        )
    print("Alignment and Tree parameters estimated.")
    print("---------------------------------------------------")

    print("Rooting gene trees and counting duplications.")
    CallDuplications.main(
        ortho_folder_path,
        output_folder_path,
        n_threads,
    )
    print("Gene tree rooting and duplication counting complete.")
    print("---------------------------------------------------")

    print("Preparing SimulationInputs...")
    sim_dir = os.path.join(output_folder_path, "TrainingResults", "SimulationInputs")
    os.makedirs(sim_dir, exist_ok=True)

    def safe_copy(src, dst_name):
        if not os.path.exists(src):
            print(f"WARNING: missing file: {src}")
            return
        shutil.copy2(src, os.path.join(sim_dir, dst_name))

    safe_copy(
        os.path.join(output_folder_path, "TrainingResults", "Orthogroup_analysis", "rescaled_ultrametric_species_tree.nwk"),
        "species_tree.nwk",
    )
    safe_copy(
        os.path.join(output_folder_path, "TrainingResults", "Orthogroup_analysis", "simulation_parameters.txt"),
        "simulation_parameters.txt",
    )
    safe_copy(
        os.path.join(output_folder_path, "TrainingResults", "Orthogroup_analysis", "gap_position_profile_counts.tsv"),
        "gap_position_profile_counts.tsv",
    )
    print(f"SimulationInputs ready at: {sim_dir}")
    print("---------------------------------------------------")

    elapsed_time = time.time() - start_time
    minutes = int(elapsed_time // 60)
    seconds = int(elapsed_time % 60)

    print("Done.")
    print(f"---Finished! This run took {minutes} minutes {seconds} seconds.")
    print(f"---Files have landed in {output_folder_path}")
    print("---Thank you for choosing OrthoSimulator.")
    print("---OrthoSimulator: Kelly Lab (2026), bioRxiv")


def main():
    parser = argparse.ArgumentParser(
        prog='OrthoFinderTrain.py',
        description=(
            "Welcome to OrthoSimulator\n"
            "Train realistic orthogroup simulations from OrthoFinder outputs.\n\n"
            "This script accepts both the OrthoSimulator.py wrapper interface and the\n"
            "old standalone interface."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    # --- New OrthoSimulator.py wrapper interface ---
    parser.add_argument(
        '--tool-output',
        dest='tool_output',
        type=str,
        default=None,
        help="Path to the OrthoFinder results folder (wrapper interface).",
    )
    parser.add_argument(
        '--output',
        dest='output',
        type=str,
        default=None,
        help="OrthoSimulator output folder (wrapper interface).",
    )
    parser.add_argument(
        '--threads',
        dest='threads',
        type=int,
        default=None,
        help="Number of threads (wrapper interface).",
    )
    # Accepted but not used — present so argparse doesn't error when the
    # outer wrapper passes them for tools that do need them.
    parser.add_argument('--tools-proteomes', dest='tools_proteomes', default=None, help=argparse.SUPPRESS)
    parser.add_argument('--tool-tree',       dest='tool_tree',       default=None, help=argparse.SUPPRESS)

    # --- Old standalone interface ---
    parser.add_argument(
        '-f',
        dest='folder',
        type=str,
        default=None,
        help="Path to the OrthoFinder results folder (standalone interface).",
    )
    parser.add_argument(
        '-o',
        dest='old_output',
        type=str,
        default=None,
        help="Output folder (standalone interface).",
    )
    parser.add_argument(
        '-t',
        dest='old_threads',
        type=int,
        default=None,
        help="Number of threads (standalone interface).",
    )

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()

    # --- Resolve args: wrapper takes priority, fall back to standalone ---
    ortho_folder_path = args.tool_output or args.folder
    output_folder_path = args.output or args.old_output
    n_threads = args.threads or args.old_threads or 8

    missing = []
    if ortho_folder_path is None:
        missing.append("--tool-output / -f (OrthoFinder results folder)")
    if output_folder_path is None:
        missing.append("--output / -o (output folder)")
    if missing:
        print("ERROR: missing required arguments:")
        for m in missing:
            print(f"  {m}")
        sys.exit(1)

    ortho_folder_path  = os.path.abspath(ortho_folder_path)
    output_folder_path = os.path.abspath(output_folder_path)

    if not os.path.isdir(ortho_folder_path):
        print(f"ERROR: OrthoFinder folder not found: {ortho_folder_path}")
        sys.exit(1)

    run(ortho_folder_path, output_folder_path, n_threads)


if __name__ == "__main__":
    main()
