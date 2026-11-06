#!/usr/bin/env python3
'''
Recompute the headline benchmark numbers from the long-format result tables.

Each table has one row per (dataset, trained_tool, tested_tool, benchmark_type,
metric). "trained_tool" is the tool whose orthogroups were used to calibrate
the simulation, and "tested_tool" is the tool that was scored on it.

For every tested tool the script prints the mean over the four calibrations and
the standard error of that mean (SD / sqrt(4)). With --ils it also compares an
ILS run against the baseline run that used the same calibration, because
comparing one ILS run with the mean of four baselines mixes two different
things.

Usage:
  python3 check_paper_numbers.py --baseline all_benchmarks_long.csv
  python3 check_paper_numbers.py --baseline all_benchmarks_long.csv \
      --ils all_benchmarks_ils_theta10.csv --ils-calibration OF3 --dataset QFO
'''

import argparse
import csv
import math
import statistics
from collections import defaultdict

HEADLINE = [
    ("orthogroup", "Recall"),
    ("orthogroup", "Precision"),
    ("orthogroup", "F1"),
    ("orthogroup", "Fission"),
    ("orthogroup", "Fusion"),
    ("ortholog", "F1"),
    ("hidden_paralogy", "rejection"),
]


def read_long(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f))


def value_table(rows, dataset):
    '''{(trained, tested, benchmark_type, metric): value}, with the derived
    hidden-paralogy rejection rate = 1 - ff_called_orthologs / ff_same_og.'''
    table = {}
    for r in rows:
        if r["dataset"] != dataset:
            continue
        table[(r["trained_tool"], r["tested_tool"], r["benchmark_type"], r["metric"])] = float(r["value"])
    derived = {}
    for (trained, tested, btype, metric), v in table.items():
        if btype == "hidden_paralogy" and metric == "ff_same_og":
            called = table.get((trained, tested, btype, "ff_called_orthologs"))
            if called is not None and v:
                derived[(trained, tested, btype, "rejection")] = 1 - called / v
    table.update(derived)
    return table


def mean_and_sem(values):
    mean = statistics.mean(values)
    sem = statistics.stdev(values) / math.sqrt(len(values)) if len(values) > 1 else float("nan")
    return mean, sem


def fmt(btype, metric, x):
    return f"{x:.2f}%" if metric in ("Fission", "Fusion") else f"{x:.4f}"


def report_baseline(table):
    tested_tools = sorted({k[1] for k in table})
    trained_tools = sorted({k[0] for k in table})
    print(f"Calibrations (trained tools): {', '.join(trained_tools)}")
    print(f"{'tested tool':<12}" + "".join(f"{b[:5]}.{m[:9]:<10}" for b, m in HEADLINE))
    for tested in tested_tools:
        cells = []
        for btype, metric in HEADLINE:
            vals = [table[(t, tested, btype, metric)] for t in trained_tools
                    if (t, tested, btype, metric) in table]
            if not vals:
                cells.append(f"{'-':<15}")
                continue
            m, s = mean_and_sem(vals)
            cells.append(f"{fmt(btype, metric, m)}±{s:.3f}".ljust(15))
        print(f"{tested:<12}" + "".join(cells))


def report_ils(base, ils, calibration, tested):
    print(f"\nILS run vs baseline with the same calibration ({calibration}), tested tool {tested}")
    print(f"{'quantity':<28}{'baseline':>12}{'ILS':>12}")
    for btype, metric in HEADLINE:
        b = base.get((calibration, tested, btype, metric))
        i = ils.get((calibration, tested, btype, metric))
        if b is None or i is None:
            continue
        print(f"{btype + ' ' + metric:<28}{fmt(btype, metric, b):>12}{fmt(btype, metric, i):>12}")
    # mean ortholog F1 / hidden-paralogy rejection over all tested tools
    for btype, metric in (("ortholog", "F1"), ("hidden_paralogy", "rejection")):
        tools = sorted({k[1] for k in base})
        bv = [base[(calibration, t, btype, metric)] for t in tools if (calibration, t, btype, metric) in base]
        iv = [ils[(calibration, t, btype, metric)] for t in tools if (calibration, t, btype, metric) in ils]
        if bv and iv:
            print(f"{'mean over tools: ' + btype + ' ' + metric:<28}"
                  f"{fmt(btype, metric, statistics.mean(bv)):>12}{fmt(btype, metric, statistics.mean(iv)):>12}")


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--baseline", required=True, help="long-format CSV of the baseline benchmarks")
    p.add_argument("--dataset", default="QFO", help="OB or QFO (default QFO)")
    p.add_argument("--ils", help="long-format CSV of an ILS run")
    p.add_argument("--ils-calibration", default="OF3", help="trained_tool used for the ILS run (default OF3)")
    p.add_argument("--ils-tested", default="of3", help="tested_tool for the single-tool ILS comparison (default of3)")
    args = p.parse_args()

    base = value_table(read_long(args.baseline), args.dataset)
    print(f"Dataset: {args.dataset}  (mean ± standard error over calibrations)")
    report_baseline(base)
    if args.ils:
        ils = value_table(read_long(args.ils), args.dataset)
        report_ils(base, ils, args.ils_calibration, args.ils_tested)


if __name__ == "__main__":
    main()
