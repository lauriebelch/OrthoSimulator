#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
FastOMA GetParameters wrapper for OrthoSimulator.

This script can be called in two ways.

1. From the outer OrthoSimulator wrapper:

    python OrthoSimulator.py \
        --Get-Parameters \
        --output Training_FastOMA \
        --tool FASTOMA \
        --tool-output /path/to/fastoma_output \
        --tools-proteomes /path/to/fastoma_input \
        --tool-tree /path/to/species_tree.nwk \
        --threads 8

2. Directly, using the old interface:

    python FastOMATrain.py \
        -fo /path/to/fastoma_output \
        -fi /path/to/fastoma_input \
        -o Training_FastOMA \
        -t 8
"""

import argparse
import os
import shutil
import sys
import time


# Ensure the script's directory is in PYTHONPATH.
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
if SCRIPT_DIR not in sys.path:
    sys.path.insert(0, SCRIPT_DIR)


import FastOMA_MakeOrthogroupSequenceFiles
import CallDuplications
import EstimateParameters


def print_train_art():
    train_art = r"""
            (  ) (@@) ( )  (@)  ()    @@    O     @     O     @
       (@@@)
   (    )
(@@@@)

     ====        ________                ___________
 _D _|  |_______/        \__I_I_____===__|_________|
  |(_)---  |   H\________/ |   |        =|___ ___|
  /     |  |   H  |  |     |   |         ||_| |_||
 |      |  |   H  |__--------------------| [___] |
 | ________|___H__/__|_____/[][]~\_______|       |
 |/ |   |-----------I_____I [][] []  D   |=======|
__/ =| o |=-~~\  /~~\  /~~\  /~~\ ____Y___________|
 |/-=|___|=O=====O=====O=====O   |_____/~\___/
  \_/      \__/  \__/  \__/  \__/

        All aboard the OrthoSimulator!
"""
    print(train_art)


def print_greeting():
    current_hour = time.localtime().tm_hour

    if 6 <= current_hour < 12:
        print("Good morning!")
    elif 12 <= current_hour < 18:
        print("Good afternoon!")
    else:
        print("Good evening!")


def safe_copy(src, dst_dir, dst_name):
    if not os.path.exists(src):
        print(f"WARNING: missing file: {src}")
        return

    os.makedirs(dst_dir, exist_ok=True)

    dst = os.path.join(dst_dir, dst_name)
    shutil.copy2(src, dst)


def validate_inputs(fastoma_out, fastoma_in, output_folder_path, n_threads):
    if not os.path.isdir(fastoma_out):
        print(f"ERROR: FastOMA output folder not found: {fastoma_out}")
        sys.exit(1)

    if not os.path.isdir(fastoma_in):
        print(f"ERROR: FastOMA input folder not found: {fastoma_in}")
        sys.exit(1)

    try:
        n_threads = int(n_threads)
    except ValueError:
        print(f"ERROR: threads must be an integer. Got: {n_threads}")
        sys.exit(1)

    if n_threads < 1:
        print(f"ERROR: threads must be at least 1. Got: {n_threads}")
        sys.exit(1)

    os.makedirs(output_folder_path, exist_ok=True)

    return n_threads


def prepare_simulation_inputs(output_folder_path):
    print("Preparing SimulationInputs...")

    sim_dir = os.path.join(
        output_folder_path,
        "TrainingResults",
        "SimulationInputs",
    )
    os.makedirs(sim_dir, exist_ok=True)

    analysis_dir = os.path.join(
        output_folder_path,
        "TrainingResults",
        "Orthogroup_analysis",
    )

    species_tree = os.path.join(
        analysis_dir,
        "rescaled_ultrametric_species_tree.nwk",
    )

    sim_params = os.path.join(
        analysis_dir,
        "simulation_parameters.txt",
    )

    gap_profile = os.path.join(
        analysis_dir,
        "gap_position_profile_counts.tsv",
    )

    safe_copy(species_tree, sim_dir, "species_tree.nwk")
    safe_copy(sim_params, sim_dir, "simulation_parameters.txt")
    safe_copy(gap_profile, sim_dir, "gap_position_profile_counts.tsv")

    print(f"SimulationInputs ready at: {sim_dir}")
    print("---------------------------------------------------")


def run_fastoma_train(
    fastoma_out,
    fastoma_in,
    output_folder_path,
    n_threads,
    tool_tree=None,
):
    fastoma_out = os.path.abspath(fastoma_out)
    fastoma_in = os.path.abspath(fastoma_in)
    output_folder_path = os.path.abspath(output_folder_path)

    if tool_tree is not None:
        tool_tree = os.path.abspath(tool_tree)

    n_threads = validate_inputs(
        fastoma_out,
        fastoma_in,
        output_folder_path,
        n_threads,
    )

    print_train_art()

    print("---------------------------------------------------")
    print("Welcome to OrthoSimulator")
    print("")
    print("Running GetParameters for FastOMA")
    print(f"FastOMA output:    {fastoma_out}")
    print(f"FastOMA input:     {fastoma_in}")
    print(f"OrthoSimulator output: {output_folder_path}")
    print(f"Threads:           {n_threads}")

    if tool_tree is not None:
        print(f"Tool tree:         {tool_tree}")

    print("")
    print_greeting()
    print("---------------------------------------------------")

    start_time = time.time()

    print("Making FastOMA sequence and alignment files ...")
    FastOMA_MakeOrthogroupSequenceFiles.main(
        fastoma_out,
        fastoma_in,
        output_folder_path,
        n_threads,
        tool_tree
    )
    print("Orthogroup sequence files created.")
    print("---------------------------------------------------")

    print("Aligning and Treeing to Estimate parameters.")
    EstimateParameters.main(
        fastoma_out,
        output_folder_path,
        n_threads,
    )
    print("Parameter estimation complete.")
    print("Alignment and Tree parameters estimated.")
    print("---------------------------------------------------")

    print("Rooting gene trees and counting duplications.")
    CallDuplications.main(
        fastoma_out,
        output_folder_path,
        n_threads,
    )
    print("Gene tree rooting and duplication counting complete.")
    print("---------------------------------------------------")

    prepare_simulation_inputs(output_folder_path)

    elapsed_time = time.time() - start_time
    minutes = int(elapsed_time // 60)
    seconds = int(elapsed_time % 60)

    print("Done.")
    print(f"---Finished! This run took {minutes} minutes {seconds} seconds.")
    print(f"---Files have landed in {output_folder_path}")
    print("---Thank you for choosing OrthoSimulator.")
    print("---OrthoSimulator: Kelly Lab (2026), bioRxiv ")


def run(tool_output, output_folder_path, n_threads, tools_proteomes=None, tool_tree=None):
    """Uniform entry point used by OrthoSimulator.py's importlib-based dispatch."""
    run_fastoma_train(tool_output, tools_proteomes, output_folder_path, n_threads, tool_tree=tool_tree)


def main():
    parser = argparse.ArgumentParser(
        prog="FastOMATrain.py",
        description=(
            "Welcome to OrthoSimulator\n"
            "Train realistic orthogroup simulations from FastOMA outputs.\n\n"
            "This script accepts both the new OrthoSimulator.py interface and the old "
            "standalone FastOMATrain.py interface."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    # New outer-wrapper interface.
    parser.add_argument(
        "--output",
        dest="output",
        type=str,
        default=None,
        help="OrthoSimulator output folder.",
    )

    parser.add_argument(
        "--tool-output",
        dest="tool_output",
        type=str,
        default=None,
        help="Path to the base FastOMA results/output folder.",
    )

    parser.add_argument(
        "--tools-proteomes",
        dest="tools_proteomes",
        type=str,
        default=None,
        help="Path to the base FastOMA input/proteomes folder.",
    )

    parser.add_argument(
        "--tool-tree",
        dest="tool_tree",
        type=str,
        default=None,
        help="Optional species tree passed by OrthoSimulator.py.",
    )

    parser.add_argument(
        "--threads",
        dest="threads",
        type=int,
        default=None,
        help="Number of threads.",
    )

    # Old standalone interface.
    parser.add_argument(
        "-fo",
        dest="fastoma_output",
        type=str,
        default=None,
        help="Path to the base FastOMA results folder.",
    )

    parser.add_argument(
        "-fi",
        dest="fastoma_input",
        type=str,
        default=None,
        help="Path to the base FastOMA input folder.",
    )

    parser.add_argument(
        "-o",
        dest="old_output",
        type=str,
        default=None,
        help="Output folder.",
    )

    parser.add_argument(
        "-t",
        dest="old_threads",
        type=int,
        default=None,
        help="Number of threads.",
    )

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()

    # Prefer the new OrthoSimulator.py arguments.
    fastoma_out = args.tool_output
    fastoma_in = args.tools_proteomes
    output_folder_path = args.output
    n_threads = args.threads

    # Fall back to old standalone arguments.
    if fastoma_out is None:
        fastoma_out = args.fastoma_output

    if fastoma_in is None:
        fastoma_in = args.fastoma_input

    if output_folder_path is None:
        output_folder_path = args.old_output

    if n_threads is None:
        n_threads = args.old_threads

    if n_threads is None:
        n_threads = 8

    missing = []

    if fastoma_out is None:
        missing.append("--tool-output or -fo")

    if fastoma_in is None:
        missing.append("--tools-proteomes or -fi")

    if output_folder_path is None:
        missing.append("--output or -o")

    if len(missing) != 0:
        print("ERROR: missing required arguments:")
        print(", ".join(missing))
        print("")
        parser.print_help()
        sys.exit(1)

    run_fastoma_train(
        fastoma_out=fastoma_out,
        fastoma_in=fastoma_in,
        output_folder_path=output_folder_path,
        n_threads=n_threads,
        tool_tree=args.tool_tree,
    )


if __name__ == "__main__":
    main()