#!/usr/bin/env python3
'''
Helpers for tests/run_long_checks.sh.

  make-config      write a simulation config that repeats a simulation from the paper
  compare-stats    compare the statistics tables of two simulations, column by column
  benchmark-values check a benchmark CSV against expected values
'''

import argparse
import csv
import os
import statistics
import sys

# columns compared between two simulations
STAT_COLUMNS = ["num_species", "num_genes", "num_duplications", "num_losses", "seq_length",
                "total_branch_length", "proportion_invariant_sites_model1", "gamma_shape_model1"]


def read_table(path):
    with open(path) as f:
        return list(csv.DictReader(f, delimiter="\t"))


def numbers(rows, column):
    out = []
    for r in rows:
        v = r.get(column)
        if v not in (None, "", "NA"):
            out.append(float(v))
    return out


##########################################
######### make-config ####################
##########################################

def cmd_make_config(args):
    '''Take the settings of a finished paper simulation (its config.txt) and the
    species tree and gap file from the matching training folder.'''
    inputs = os.path.join(args.training, "TrainingResults", "SimulationInputs")
    skip = {"threads", "SimPath", "Orthogroups", "ultrametric_tree", "gap_file"}
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    with open(args.paper_config) as f, open(args.out, "w") as w:
        for line in f:
            if "=" in line and line.split("=", 1)[0].strip() not in skip:
                w.write(line.rstrip("\n") + "\n")
        w.write(f"ultrametric_tree={os.path.join(inputs, 'species_tree.nwk')}\n")
        w.write(f"gap_file={os.path.join(inputs, 'gap_position_profile_counts.tsv')}\n")
        w.write(f"Orthogroups={args.n}\n")
    print(f"wrote {args.out}")


##########################################
######### compare-stats ##################
##########################################

def cmd_compare_stats(args):
    from scipy import stats
    a, b = read_table(args.a), read_table(args.b)
    print(f"A = {args.a} ({len(a)} orthogroups)\nB = {args.b} ({len(b)} orthogroups)")
    print(f"{'column':<36}{'mean A':>11}{'mean B':>11}{'median A':>10}{'median B':>10}{'KS D':>8}{'p':>8}")
    smallest_p = 1.0
    for col in STAT_COLUMNS:
        x, y = numbers(a, col), numbers(b, col)
        if not x or not y:
            continue
        ks = stats.ks_2samp(x, y)
        smallest_p = min(smallest_p, ks.pvalue)
        print(f"{col:<36}{statistics.mean(x):11.3f}{statistics.mean(y):11.3f}"
              f"{statistics.median(x):10.2f}{statistics.median(y):10.2f}{ks.statistic:8.3f}{ks.pvalue:8.3f}")
    print(f"smallest p = {smallest_p:.2g}   (fails below {args.min_p})")
    if smallest_p < args.min_p:
        print("FAIL  a distribution differs more than chance would explain")
        sys.exit(1)
    print("PASS")


##########################################
######### benchmark-values ###############
##########################################

def cmd_benchmark_values(args):
    with open(args.csv) as f:
        row = next(csv.DictReader(f))
    failed = 0
    for item in args.expect:
        name, value = item.split("=")
        got = float(row[name])
        ok = abs(got - float(value)) <= args.tolerance
        print(f"  {'PASS' if ok else 'FAIL'}  {name} = {got} (expected {value} +/- {args.tolerance})")
        failed += not ok
    sys.exit(1 if failed else 0)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="command", required=True)

    q = sub.add_parser("make-config")
    q.add_argument("--paper-config", required=True, help="config.txt of a simulation from the paper")
    q.add_argument("--training", required=True, help="the matching Training_<tool>_<dataset> folder")
    q.add_argument("--n", type=int, required=True)
    q.add_argument("--out", required=True)
    q.set_defaults(func=cmd_make_config)

    q = sub.add_parser("compare-stats")
    q.add_argument("--a", required=True)
    q.add_argument("--b", required=True)
    q.add_argument("--min-p", type=float, default=1e-4,
                   help="fail if any column's KS p-value is below this (default 1e-4)")
    q.set_defaults(func=cmd_compare_stats)

    q = sub.add_parser("benchmark-values")
    q.add_argument("--csv", required=True)
    q.add_argument("--expect", nargs="+", required=True, help="Column=value pairs")
    q.add_argument("--tolerance", type=float, default=1e-4)
    q.set_defaults(func=cmd_benchmark_values)

    args = p.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
