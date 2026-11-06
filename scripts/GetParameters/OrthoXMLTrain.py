#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
OrthoXML GetParameters wrapper for OrthoSimulator.

Generic input adapter for any orthology tool that can export (or be
converted to) the OrthoXML standard, rather than a bespoke wrapper per
tool. See OrthoXML_MakeOrthogroupSequenceFiles.py for the scope decision:
OrthoXML is only used for root-level orthogroup membership - alignment,
tree estimation, and duplication calling all run identically to every
other tool in this pipeline.

Usage (via the outer OrthoSimulator.py wrapper):

    python OrthoSimulator.py \
        --Get-Parameters \
        --output Training_OrthoXML \
        --tool ORTHOXML \
        --tool-output /path/to/results.orthoxml \
        --tools-proteomes /path/to/proteomes_dir \
        --tool-tree /path/to/species_tree.nwk \
        --threads 8
"""

import os
import sys
import time


# Ensure the script's directory is in PYTHONPATH before any local imports.
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
if SCRIPT_DIR not in sys.path:
    sys.path.insert(0, SCRIPT_DIR)

import OrthoXML_MakeOrthogroupSequenceFiles
import CallDuplications
import EstimateParameters


def print_greeting():
    current_hour = time.localtime().tm_hour
    if 6 <= current_hour < 12:
        print("Good morning!")
    elif 12 <= current_hour < 18:
        print("Good afternoon!")
    else:
        print("Good evening!")


def validate_inputs(orthoxml_file, proteome_dir, tool_tree, output_folder_path, n_threads):
    if not os.path.isfile(orthoxml_file):
        print(f"ERROR: OrthoXML file not found: {orthoxml_file}")
        sys.exit(1)

    if not proteome_dir or not os.path.isdir(proteome_dir):
        print(f"ERROR: proteome directory not found: {proteome_dir}")
        sys.exit(1)

    if not tool_tree or not os.path.isfile(tool_tree):
        print(f"ERROR: species tree not found: {tool_tree}")
        sys.exit(1)

    try:
        n_threads = int(n_threads)
    except (TypeError, ValueError):
        print(f"ERROR: threads must be an integer. Got: {n_threads}")
        sys.exit(1)

    if n_threads < 1:
        print(f"ERROR: threads must be at least 1. Got: {n_threads}")
        sys.exit(1)

    os.makedirs(output_folder_path, exist_ok=True)

    return n_threads


def run_orthoxml_train(
    orthoxml_file,
    tools_proteomes,
    output_folder_path,
    n_threads,
    tool_tree=None,
):
    orthoxml_file = os.path.abspath(orthoxml_file)
    proteome_dir = os.path.abspath(tools_proteomes) if tools_proteomes else None
    output_folder_path = os.path.abspath(output_folder_path)
    tool_tree_abs = os.path.abspath(tool_tree) if tool_tree else None

    n_threads = validate_inputs(
        orthoxml_file, proteome_dir, tool_tree_abs, output_folder_path, n_threads
    )

    print("---------------------------------------------------")
    print("Welcome to OrthoSimulator")
    print("")
    print("Running GetParameters for OrthoXML")
    print(f"OrthoXML file:     {orthoxml_file}")
    print(f"Proteomes:         {proteome_dir}")
    print(f"Species tree:      {tool_tree_abs}")
    print(f"OrthoSimulator output: {output_folder_path}")
    print(f"Threads:           {n_threads}")
    print("")
    print_greeting()
    print("---------------------------------------------------")

    start_time = time.time()

    print("Making OrthoXML-derived orthogroup sequence files ...")
    OrthoXML_MakeOrthogroupSequenceFiles.main(
        orthoxml_file,
        proteome_dir,
        output_folder_path,
        n_threads,
        tool_tree_abs,
    )
    print("Orthogroup sequence files created.")
    print("---------------------------------------------------")

    print("Aligning and Treeing to Estimate parameters.")
    EstimateParameters.main(
        output_folder_path,
        output_folder_path,
        n_threads,
    )
    print("Parameter estimation complete.")
    print("Alignment and Tree parameters estimated.")
    print("---------------------------------------------------")

    print("Rooting gene trees and counting duplications.")
    # NOTE: output_folder_path is passed as the "tool output directory" here
    # (rather than orthoxml_file, which is a single file with no directory
    # of its own) because OrthoXML_MakeOrthogroupSequenceFiles.py writes
    # Species_Tree/ under output_folder_path - see that script's main() for
    # why. CallDuplications.py just needs a directory containing
    # Species_Tree/SpeciesTree_rooted_node_labels_numeric.txt.
    CallDuplications.main(
        output_folder_path,
        output_folder_path,
        n_threads,
    )
    print("Gene tree rooting and duplication counting complete.")
    print("---------------------------------------------------")

    elapsed_time = time.time() - start_time
    minutes = int(elapsed_time // 60)
    seconds = int(elapsed_time % 60)

    print("Done.")
    print(f"---Finished! This run took {minutes} minutes {seconds} seconds.")
    print(f"---Files have landed in {output_folder_path}")
    print("---Thank you for choosing OrthoSimulator.")


def run(tool_output, output_folder_path, n_threads, tools_proteomes=None, tool_tree=None):
    """Uniform entry point used by OrthoSimulator.py's importlib-based dispatch."""
    run_orthoxml_train(
        tool_output,
        tools_proteomes,
        output_folder_path,
        n_threads,
        tool_tree=tool_tree,
    )


if __name__ == "__main__":
    print(
        "OrthoXMLTrain.py is intended to be run via OrthoSimulator.py --Get-Parameters "
        "--tool ORTHOXML. See the module docstring for the expected flags."
    )
    sys.exit(1)
