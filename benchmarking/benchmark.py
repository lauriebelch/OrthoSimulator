#!/usr/bin/env python3
"""
benchmark.py — Orthogroup benchmarking against a simulated reference.

Usage:
    python benchmark.py --output-folder <sim_dir> --tool-output <tool_dir> --tool <tool>
                        [--csv <results.csv>]

Arguments:
    --output-folder / -o    Simulation output folder (contains Simulated_Orthogroups.txt).
    --tool-output   / -to   The tool's output folder.
    --tool          / -t    Tool name: orthofinder, fastoma, broccoli, sonicparanoid2
    --csv                   Optional path to write results as a CSV file.

Tool-specific file locations (relative to --tool-output):
    orthofinder   : **/Orthogroups/Orthogroups.txt  (searched recursively)
    fastoma       : RootHOGs.tsv
    broccoli      : dir_step3/orthologous_groups.txt
    sonicparanoid2: sp_default/runs/*/ortholog_groups/flat.ortholog_groups.tsv
"""

import argparse
import csv
import glob
import math
import os
import sys
from collections import defaultdict
from pathlib import Path


SUPPORTED_TOOLS = ("orthofinder", "fastoma", "broccoli", "sonicparanoid2")


# ---------------------------------------------------------------------------
# Tool-specific loaders
# ---------------------------------------------------------------------------

def find_orthofinder_og_file(tool_output):
    matches = glob.glob(
        str(tool_output / "**" / "Orthogroups" / "Orthogroups.txt"),
        recursive=True,
    )
    if not matches:
        raise FileNotFoundError(
            f"Could not find Orthogroups/Orthogroups.txt under {tool_output}.\n"
            "Expected path: <tool_output>/<any>/Results_*/Orthogroups/Orthogroups.txt"
        )
    if len(matches) > 1:
        matches = sorted(matches, key=os.path.getmtime, reverse=True)
        print(f"Multiple Orthogroups.txt found; using most recent: {matches[0]}")
    return Path(matches[0])


def load_orthofinder(file_path):
    og_tbl = []
    with open(file_path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            og_id = line.split(":")[0]
            genes = line.split(":")[1].strip().split()
            for g in genes:
                gene_clean = g.split(".")[0]
                fields = gene_clean.split("_")
                if len(fields) < 3:
                    print(f"Skipping unparseable gene: {gene_clean}")
                    continue
                try:
                    refog = int(fields[2])
                except ValueError:
                    print(f"Skipping gene with non-numeric refOG field: {gene_clean}")
                    continue
                og_tbl.append((og_id, gene_clean, refog))
    return og_tbl


def load_fastoma(file_path):
    og_tbl = []
    with open(file_path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            og_id = parts[0].rstrip(":").strip()
            gene = parts[1].strip()
            if gene.lower() in {"protein", "gene", "genes"}:
                continue
            if og_id.lower() in {"hog", "roothog", "orthogroup"}:
                continue
            if not gene:
                continue
            gene_clean = gene.split(".")[0]
            fields = gene_clean.split("_")
            if len(fields) < 3:
                print(f"Skipping unparseable gene row: {line}")
                continue
            try:
                refog = int(fields[2])
            except ValueError:
                print(f"Skipping gene with non-numeric refOG field: {gene_clean}")
                continue
            og_tbl.append((og_id, gene_clean, refog))
    return og_tbl


def load_broccoli(file_path):
    og_tbl = []
    with open(file_path) as f:
        next(f, None)  # skip header
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            og_id = parts[0]
            genes = parts[1:]
            for gene in genes:
                if not gene or gene == "*":
                    continue
                gene_clean = gene.split(".")[0]
                fields = gene_clean.split("_")
                if len(fields) < 3:
                    print(f"Skipping unparseable gene: {gene_clean}")
                    continue
                try:
                    refog = int(fields[2])
                except ValueError:
                    print(f"Skipping gene with non-numeric refOG field: {gene_clean}")
                    continue
                og_tbl.append((og_id, gene_clean, refog))
    return og_tbl


def find_sonicparanoid_og_file(tool_output):
    for base in [tool_output / "sp_default", tool_output]:
        pattern = str(base / "runs" / "*" / "ortholog_groups" / "flat.ortholog_groups.tsv")
        matches = glob.glob(pattern)
        if matches:
            break
    if not matches:
        raise FileNotFoundError(
            f"No flat.ortholog_groups.tsv found under {tool_output}.\n"
            "Expected: <tool_output>/sp_default/runs/*/ortholog_groups/flat.ortholog_groups.tsv"
        )
    if len(matches) > 1:
        matches = sorted(matches, key=os.path.getmtime, reverse=True)
        print(f"Multiple SonicParanoid2 runs found; using most recent: {matches[0]}")
    return Path(matches[0])


def load_sonicparanoid(file_path):
    og_tbl = []
    with open(file_path) as f:
        header = f.readline().rstrip("\n").split("\t")
        if "group_id" not in header:
            raise ValueError(f"'group_id' column not found in: {file_path}")
        group_idx = header.index("group_id")
        species_indices = [i for i, col in enumerate(header) if i != group_idx]
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < len(header):
                parts = parts + [""] * (len(header) - len(parts))
            og_id = parts[group_idx].strip()
            if not og_id:
                continue
            for idx in species_indices:
                cell = parts[idx].strip()
                if not cell or cell == "*" or cell.upper() == "NA":
                    continue
                genes = [g.strip() for g in cell.split(",") if g.strip() and g.strip() != "*"]
                for gene in genes:
                    gene_clean = gene.split(".")[0]
                    fields = gene_clean.split("_")
                    if len(fields) < 3:
                        print(f"Skipping unparseable gene: {gene_clean}")
                        continue
                    try:
                        refog = int(fields[2])
                    except ValueError:
                        print(f"Skipping gene with non-numeric refOG field: {gene_clean}")
                        continue
                    og_tbl.append((og_id, gene_clean, refog))
    return og_tbl


# ---------------------------------------------------------------------------
# Shared metrics
# ---------------------------------------------------------------------------

def remove_singletons(og_tbl):
    counts = defaultdict(int)
    for og, _, _ in og_tbl:
        counts[og] += 1
    return [row for row in og_tbl if counts[row[0]] > 1]


def load_reference(file_path):
    ref_genes = []
    ref_og_ids = set()
    with open(file_path) as f:
        f.readline()  # skip header
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            og = parts[0]
            ref_og_ids.add(int(og))
            for cell in parts[1:]:
                if cell:
                    for g in cell.split(","):
                        g = g.strip()
                        if g:
                            ref_genes.append((int(og), g))
    return ref_og_ids, ref_genes


def compute_missing(ref_og_ids, ref_genes, og_tbl):
    pred_refogs = {r for _, _, r in og_tbl}
    pred_genes  = {g for _, g, _ in og_tbl}
    missing_refog = len([r for r in ref_og_ids if r not in pred_refogs]) / len(ref_og_ids) * 100
    missing_genes = len([g for _, g in ref_genes if g not in pred_genes]) / len(ref_genes) * 100
    return missing_refog, missing_genes


def compute_fusion(ref_og_ids, og_tbl):
    fused = set()
    og_map = defaultdict(list)
    for og, gene, refog in og_tbl:
        og_map[og].append(refog)
    for og in og_map:
        if len(set(og_map[og])) > 1:
            fused.update(og_map[og])
    return len(fused) / len(ref_og_ids) * 100


def compute_fission(ref_og_ids, og_tbl):
    ref_map = defaultdict(list)
    for og, gene, refog in og_tbl:
        ref_map[refog].append(og)
    fissed = [r for r in ref_og_ids if len(set(ref_map[r])) > 1]
    return len(fissed) / len(ref_og_ids) * 100


def compute_recall_precision(ref_og_ids, ref_genes, og_tbl):
    ref_map  = defaultdict(set)
    for og, gene in ref_genes:
        ref_map[og].add(gene)
    pred_map = defaultdict(set)
    for og, gene, _ in og_tbl:
        pred_map[og].add(gene)

    RC_num_rec = RC_den_rec = 0
    RC_num_pre = RC_den_pre = 0

    for refog in ref_og_ids:
        ref_set   = ref_map[refog]
        size_Ri   = len(ref_set)
        num_i_rec = den_i_rec = 0
        num_i_pre = den_i_pre = 0
        for predog in pred_map:
            pred_set = pred_map[predog]
            overlap  = len(ref_set & pred_set)
            if overlap == 0:
                continue
            rc_rec = overlap / size_Ri
            num_i_rec += overlap * rc_rec
            den_i_rec += overlap
            rc_pre = overlap / len(pred_set)
            num_i_pre += overlap * rc_pre
            den_i_pre += overlap
        RC_i_rec = num_i_rec / den_i_rec if den_i_rec else 0
        RC_i_pre = num_i_pre / den_i_pre if den_i_pre else 0
        RC_num_rec += size_Ri * RC_i_rec
        RC_den_rec += size_Ri
        RC_num_pre += size_Ri * RC_i_pre
        RC_den_pre += size_Ri

    recall    = RC_num_rec / RC_den_rec if RC_den_rec else 0
    precision = RC_num_pre / RC_den_pre if RC_den_pre else 0
    return recall, precision


def compute_entropy(ref_og_ids, ref_genes, og_tbl):
    ref_map = defaultdict(list)
    for og, gene in ref_genes:
        ref_map[og].append(gene)

    total_entropy = 0
    total_genes   = 0

    for refog in ref_og_ids:
        genes  = ref_map[refog]
        size   = len(genes)
        counts = defaultdict(int)
        for predog, gene, r in og_tbl:
            if r == refog:
                counts[predog] += 1
        entropy_val = 0
        for c in counts.values():
            p = c / size
            entropy_val += -p * math.log2(p)
        if size > 1:
            entropy_val /= math.log2(size)
        total_entropy += entropy_val * size
        total_genes   += size

    return total_entropy / total_genes if total_genes else 0


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description=(
            "Benchmark orthogroup predictions against a simulated reference.\n\n"
            "Requires two folders:\n"
            "  --output-folder  : simulation folder containing Simulated_Orthogroups.txt\n"
            "  --tool-output    : the tool's output folder"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--output-folder", "-o", required=True, dest="output_folder",
                        help="Simulation output folder (contains Simulated_Orthogroups.txt).")
    parser.add_argument("--tool-output", "-to", required=True, dest="tool_output",
                        help="Tool output folder.")
    parser.add_argument("--tool", "-t", required=True, dest="tool", choices=SUPPORTED_TOOLS,
                        help=f"Orthology tool. Supported: {', '.join(SUPPORTED_TOOLS)}")
    parser.add_argument("--csv", default=None, dest="csv",
                        help="Optional path to write results as a CSV file.")

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()

    output_folder = Path(args.output_folder).resolve()
    tool_output   = Path(args.tool_output).resolve()
    tool          = args.tool.lower()

    for label, folder in [("output-folder", output_folder), ("tool-output", tool_output)]:
        if not folder.is_dir():
            print(f"ERROR: --{label} not found: {folder}")
            sys.exit(1)

    ref_file = output_folder / "Simulated_Orthogroups.txt"
    if not ref_file.exists():
        print(f"ERROR: reference file not found: {ref_file}")
        sys.exit(1)

    print(f"Tool:            {tool}")
    print(f"Simulation dir:  {output_folder}")
    print(f"Tool output dir: {tool_output}")

    # --- Load tool predictions ---
    if tool == "orthofinder":
        og_file = find_orthofinder_og_file(tool_output)
        print(f"OG file:         {og_file}")
        og_tbl = load_orthofinder(og_file)
    elif tool == "fastoma":
        og_file = tool_output / "RootHOGs.tsv"
        if not og_file.exists():
            print(f"ERROR: FastOMA file not found: {og_file}"); sys.exit(1)
        print(f"OG file:         {og_file}")
        og_tbl = load_fastoma(og_file)
    elif tool == "broccoli":
        og_file = tool_output / "dir_step3" / "orthologous_groups.txt"
        if not og_file.exists():
            print(f"ERROR: Broccoli file not found: {og_file}"); sys.exit(1)
        print(f"OG file:         {og_file}")
        og_tbl = load_broccoli(og_file)
    elif tool == "sonicparanoid2":
        og_file = find_sonicparanoid_og_file(tool_output)
        print(f"OG file:         {og_file}")
        og_tbl = load_sonicparanoid(og_file)

    og_tbl = remove_singletons(og_tbl)
    ref_og_ids, ref_genes = load_reference(ref_file)

    # --- Compute metrics ---
    missing_refog, missing_genes = compute_missing(ref_og_ids, ref_genes, og_tbl)
    fusion    = compute_fusion(ref_og_ids, og_tbl)
    fission   = compute_fission(ref_og_ids, og_tbl)
    recall, precision = compute_recall_precision(ref_og_ids, ref_genes, og_tbl)
    f1        = 2 * recall * precision / (recall + precision) if (recall + precision) > 0 else 0
    entropy   = compute_entropy(ref_og_ids, ref_genes, og_tbl)

    results = {
        "tool":             tool,
        "simulation_folder": str(output_folder),
        "tool_output":      str(tool_output),
        "Missing_Refog":    round(missing_refog, 6),
        "Missing_Genes":    round(missing_genes, 6),
        "Fusion":           round(fusion, 6),
        "Fission":          round(fission, 6),
        "Recall":           round(recall, 6),
        "Precision":        round(precision, 6),
        "F1":               round(f1, 6),
        "Entropy":          round(entropy, 6),
    }

    # --- Print results ---
    print("\nRESULTS\n")
    print(f"Missing_Refog (%): {missing_refog:.4f}")
    print(f"Missing_Genes (%): {missing_genes:.4f}")
    print(f"Fusion        (%): {fusion:.4f}")
    print(f"Fission       (%): {fission:.4f}")
    print(f"Recall:            {recall:.4f}")
    print(f"Precision:         {precision:.4f}")
    print(f"F1:                {f1:.4f}")
    print(f"Entropy:           {entropy:.4f}")

    # --- Write CSV ---
    if args.csv:
        csv_path = Path(args.csv)
        csv_path.parent.mkdir(parents=True, exist_ok=True)
        with open(csv_path, "w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=list(results.keys()))
            writer.writeheader()
            writer.writerow(results)
        print(f"\nResults written to: {csv_path}")


if __name__ == "__main__":
    main()
