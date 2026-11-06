#!/usr/bin/env python3
"""
bootstrap_rank_summarize.py — step 2: combine the four per-tool CSVs written
by bootstrap_rank_support.py and report how well-supported the published
OrthoBench tool ranking is under repeated resampling.

Run bootstrap_rank_support.py once per tool first (same --samples-file each
time), then:

    python bootstrap_rank_summarize.py \
        --orthofinder-csv bootstrap_metrics_orthofinder.csv \
        --fastoma-csv bootstrap_metrics_fastoma.csv \
        --broccoli-csv bootstrap_metrics_broccoli.csv \
        --sonicparanoid2-csv bootstrap_metrics_sonicparanoid2.csv \
        --published-tsv orthobench_published.tsv \
        --csv bootstrap_rank_support_results.csv

For each metric: what fraction of replicates reproduce a tool ranking
consistent with the published order (ties in the published data are treated
as blocks - within-block order doesn't need to match), and whether the
published order is the single most frequently observed ranking (and if not,
what rank by frequency it occupies).
"""

import argparse
import csv
import importlib.util
import sys
from collections import Counter
from pathlib import Path


def load_bootstrap_module():
    path = Path(__file__).resolve().parent / "bootstrap_rank_support.py"
    spec = importlib.util.spec_from_file_location("bootstrap_rank_support", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def read_metrics_csv(path):
    meta = {}
    rows = []
    with open(path, newline="") as f:
        lines = f.readlines()
    data_lines = []
    for line in lines:
        if line.startswith("#"):
            key, _, value = line[1:].strip().partition("=")
            meta[key.strip()] = value.strip()
        else:
            data_lines.append(line)
    reader = csv.DictReader(data_lines)
    for row in reader:
        rows.append({
            "replicate": int(row["replicate"]),
            "Missing_Refog": float(row["Missing_Refog"]),
            "Missing_Genes": float(row["Missing_Genes"]),
            "Fusion": float(row["Fusion"]),
            "Fission": float(row["Fission"]),
            "Recall": float(row["Recall"]),
            "Precision": float(row["Precision"]),
            "Entropy": float(row["Entropy"]),
        })
    return meta, rows


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--orthofinder-csv", required=True)
    parser.add_argument("--fastoma-csv", required=True)
    parser.add_argument("--broccoli-csv", required=True)
    parser.add_argument("--sonicparanoid2-csv", required=True)
    parser.add_argument("--published-tsv", required=True)
    parser.add_argument("--of3-row", default=None)
    parser.add_argument("--fastoma-row", default=None)
    parser.add_argument("--broccoli-row", default=None)
    parser.add_argument("--sp2-row", default=None)
    parser.add_argument("--csv", default=None,
                        help="Optional path to write the results table as CSV.")

    args = parser.parse_args()

    brs = load_bootstrap_module()
    TOOLS = brs.TOOLS
    METRICS = brs.METRICS
    DEFAULT_ROWS = brs.DEFAULT_ROWS

    row_map = {
        "orthofinder": args.of3_row or DEFAULT_ROWS["orthofinder"],
        "fastoma": args.fastoma_row or DEFAULT_ROWS["fastoma"],
        "broccoli": args.broccoli_row or DEFAULT_ROWS["broccoli"],
        "sonicparanoid2": args.sp2_row or DEFAULT_ROWS["sonicparanoid2"],
    }

    csv_paths = {
        "orthofinder": args.orthofinder_csv,
        "fastoma": args.fastoma_csv,
        "broccoli": args.broccoli_csv,
        "sonicparanoid2": args.sonicparanoid2_csv,
    }

    meta = {}
    rows = {}
    for tool in TOOLS:
        meta[tool], rows[tool] = read_metrics_csv(csv_paths[tool])
        if meta[tool].get("tool") != tool:
            print(f"WARNING: {csv_paths[tool]} was recorded for tool "
                  f"'{meta[tool].get('tool')}', not '{tool}'. Check you passed "
                  "the right files to the right flags.")

    # Sanity check all four tools used the same samples file and replicate count.
    samples_files = {meta[t].get("samples_file") for t in TOOLS}
    if len(samples_files) > 1:
        print("WARNING: the four CSVs were generated with different "
              f"--samples-file values ({samples_files}). Per-replicate "
              "rankings across tools are only meaningful if all four used "
              "the exact same subsamples - results below may not be valid.")

    n_replicates = {len(rows[t]) for t in TOOLS}
    if len(n_replicates) > 1:
        print(f"ERROR: tools have different replicate counts: "
              f"{ {t: len(rows[t]) for t in TOOLS} }")
        sys.exit(1)
    n = n_replicates.pop()

    for t in TOOLS:
        reps = [r["replicate"] for r in rows[t]]
        if reps != list(range(n)):
            print(f"ERROR: {csv_paths[t]} does not have a clean 0..{n-1} "
                  "replicate index sequence - cannot align tools by replicate.")
            sys.exit(1)

    published_values = brs.load_published_values(args.published_tsv, row_map)
    published_blocks = {
        metric: brs.build_published_blocks(published_values, column, direction, TOOLS)
        for metric, (column, direction) in METRICS.items()
    }

    print("Published order (best -> worst), from the table:")
    for metric, blocks in published_blocks.items():
        block_strs = [
            "{" + ", ".join(sorted(b)) + "}" if len(b) > 1 else next(iter(b))
            for b in blocks
        ]
        print(f"  {metric:15s}: " + " > ".join(block_strs))
    print("")

    match_counts = {metric: 0 for metric in METRICS}
    ordering_counts = {metric: Counter() for metric in METRICS}

    for rep in range(n):
        values_by_metric = {metric: {} for metric in METRICS}
        for metric in METRICS:
            for tool in TOOLS:
                values_by_metric[metric][tool] = rows[tool][rep][metric]

        for metric, (_, direction) in METRICS.items():
            ordering = brs.rank_order(values_by_metric[metric], direction, TOOLS)
            ordering_counts[metric][ordering] += 1
            if brs.consistent_with_blocks(ordering, published_blocks[metric]):
                match_counts[metric] += 1

    print("RESULTS")
    print("")

    results_rows = []
    for metric in METRICS:
        published_bucket_count = match_counts[metric]
        match_pct = published_bucket_count / n * 100

        entries = [("<published order>", published_bucket_count)]
        for ordering, count in ordering_counts[metric].items():
            if not brs.consistent_with_blocks(ordering, published_blocks[metric]):
                entries.append((ordering, count))
        entries.sort(key=lambda e: e[1], reverse=True)

        published_rank = 1 + sum(
            1 for label, count in entries
            if label != "<published order>" and count > published_bucket_count
        )
        is_number_one = (published_rank == 1)

        top_label, top_count = entries[0]
        top_is_published = (top_label == "<published order>")
        top_str = "<published order>" if top_is_published else " > ".join(top_label)

        n_distinct_orderings = len(ordering_counts[metric])

        print(f"{metric}:")
        print(f"  Published order reproduced in {published_bucket_count}/{n} "
              f"resamples ({match_pct:.1f}%)")
        print(f"  Rank of published order by frequency: #{published_rank} "
              f"(out of {len(entries)} distinct patterns)")
        print(f"  Most frequent pattern overall: {top_str} "
              f"({top_count}/{n} = {top_count / n * 100:.1f}%)")
        print(f"  Distinct full orderings observed: {n_distinct_orderings}")
        print("")

        results_rows.append({
            "metric": metric,
            "published_match_count": published_bucket_count,
            "replicates": n,
            "published_match_pct": round(match_pct, 2),
            "published_rank_by_frequency": published_rank,
            "is_best_supported": is_number_one,
            "n_distinct_full_orderings": n_distinct_orderings,
            "top_pattern_is_published": top_is_published,
            "top_pattern_count": top_count,
            "top_pattern_pct": round(top_count / n * 100, 2),
        })

    if args.csv:
        csv_path = Path(args.csv)
        csv_path.parent.mkdir(parents=True, exist_ok=True)
        with open(csv_path, "w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=list(results_rows[0].keys()))
            writer.writeheader()
            writer.writerows(results_rows)
        print(f"Results written to: {csv_path}")


if __name__ == "__main__":
    main()
