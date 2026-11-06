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

import Broccoli_MakeOrthogroupSequenceFiles
import CallDuplications
import EstimateParameters


def run(tool_output, output_folder_path, n_threads, tools_proteomes=None, tool_tree=None):
    if tools_proteomes is None:
        raise ValueError("BroccoliTrain.run(): tools_proteomes is required but was not provided.")
    if tool_tree is None:
        raise ValueError("BroccoliTrain.run(): tool_tree is required but was not provided.")

    broccoli_out       = os.path.abspath(tool_output)
    broccoli_in        = os.path.abspath(tools_proteomes)
    output_folder_path = os.path.abspath(output_folder_path)
    tool_tree          = os.path.abspath(tool_tree)

    os.makedirs(output_folder_path, exist_ok=True)

    # --- Welcome ---
    print("---------------------------------------------------")
    print("Welcome to OrthoSimulator\n")
    print("Running GetParameters for Broccoli")
    print(f"Broccoli output:   {broccoli_out}")
    print(f"Broccoli input:    {broccoli_in}")
    print(f"Species tree:      {tool_tree}")
    print(f"OrthoSimulator output: {output_folder_path}")
    print(f"Threads:           {n_threads}")
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
    print("Making Broccoli sequence files ...")
    Broccoli_MakeOrthogroupSequenceFiles.main(
        broccoli_out,
        broccoli_in,
        output_folder_path,
        n_threads,
        tool_tree,
    )
    print("Orthogroup sequence files created.")
    print("---------------------------------------------------")

    print("Aligning and Treeing to Estimate parameters.")
    EstimateParameters.main(
        broccoli_out,
        output_folder_path,
        n_threads,
    )
    print("Parameter estimation complete.")
    print("---------------------------------------------------")

    print("Rooting gene trees and counting duplications.")
    CallDuplications.main(
        broccoli_out,
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
        prog='BroccoliTrain.py',
        description=(
            "Welcome to OrthoSimulator\n"
            "Train realistic orthogroup simulations from Broccoli outputs.\n\n"
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
        help="Path to the Broccoli results folder (wrapper interface).",
    )
    parser.add_argument(
        '--tools-proteomes',
        dest='tools_proteomes',
        type=str,
        default=None,
        help="Path to the Broccoli input proteomes folder (wrapper interface).",
    )
    parser.add_argument(
        '--tool-tree',
        dest='tool_tree',
        type=str,
        default=None,
        help="Path to the species tree .nwk file (required — Broccoli does not produce one natively).",
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

    # --- Old standalone interface ---
    parser.add_argument(
        '-bo',
        dest='broccoli_output',
        type=str,
        default=None,
        help="Path to the Broccoli results folder (standalone interface).",
    )
    parser.add_argument(
        '-bi',
        dest='broccoli_input',
        type=str,
        default=None,
        help="Path to the Broccoli input proteomes folder (standalone interface).",
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
    broccoli_out = args.tool_output or args.broccoli_output
    broccoli_in  = args.tools_proteomes or args.broccoli_input
    output_folder_path = args.output or args.old_output
    n_threads = args.threads or args.old_threads or 8
    tool_tree = args.tool_tree

    # --- Validate required args ---
    missing = []
    if broccoli_out is None:
        missing.append("--tool-output / -bo (Broccoli results folder)")
    if broccoli_in is None:
        missing.append("--tools-proteomes / -bi (Broccoli input proteomes folder)")
    if output_folder_path is None:
        missing.append("--output / -o (output folder)")
    if tool_tree is None:
        missing.append("--tool-tree (species tree — Broccoli does not produce one natively)")
    if missing:
        print("ERROR: missing required arguments:")
        for m in missing:
            print(f"  {m}")
        sys.exit(1)

    broccoli_out = os.path.abspath(broccoli_out)
    broccoli_in  = os.path.abspath(broccoli_in)
    output_folder_path = os.path.abspath(output_folder_path)
    tool_tree = os.path.abspath(tool_tree)

    if not os.path.isdir(broccoli_out):
        print(f"ERROR: Broccoli output folder not found: {broccoli_out}")
        sys.exit(1)
    if not os.path.isdir(broccoli_in):
        print(f"ERROR: Broccoli input proteomes folder not found: {broccoli_in}")
        sys.exit(1)
    if not os.path.isfile(tool_tree):
        print(f"ERROR: Species tree not found: {tool_tree}")
        sys.exit(1)

    run(broccoli_out, output_folder_path, n_threads, tools_proteomes=broccoli_in, tool_tree=tool_tree)


if __name__ == "__main__":
    main()
