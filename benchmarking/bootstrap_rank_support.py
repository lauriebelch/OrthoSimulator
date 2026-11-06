#!/usr/bin/env python3
"""
bootstrap_rank_support.py — test how well-supported the published OrthoBench
tool ranking is, by repeated resampling from a large simulated reference.

What it does
------------
The published OrthoBench comparison (orthobench_published.tsv) ranks four
tools using only 70 reference orthogroups. Our simulations give complete
ground truth for ~10,000, so we can ask a sharper question than "do the
numbers match": if OrthoBench had happened to curate a different random 70
orthogroups from the same population, how often would the same tool
*ranking* come out? This script:

  1. Reads the published order for the four tools (OrthoFinder, FastOMA,
     Broccoli, SonicParanoid2) from orthobench_published.tsv, per metric,
     respecting ties in the published values (tools tied in the real data
     are treated as a block - within-block order doesn't need to match).
  2. Draws N random samples of --sample-size reference orthogroups (without
     replacement) from the simulation's full reference set.
  3. Recomputes all four tools' metrics on each sample using the exact same
     per-orthogroup logic as benchmark.py (imported directly, not
     duplicated), and ranks the four tools by each metric.
  4. Reports, per metric: what fraction of resamples produce a ranking
     consistent with the published order, and whether the published order
     is the single most frequently observed ranking overall (and if not,
     what rank by frequency it occupies).

Only rank/order is compared - never raw magnitude. A tool being simulated
at 99% vs a published 100% counts identically to 60% vs 61%, as long as the
relative order of all four tools is what matters.

No Kendall's tau, no equivalence thresholds - published ties are handled as
block constraints (does the resampled ranking respect the between-block
order, ignoring how tied tools sort against each other), and "is it #1" is
answered by literally counting and ranking the distinct orderings observed.

Usage:
    python bootstrap_rank_support.py \
        --output-folder <sim_dir> \
        --orthofinder-output <dir> --fastoma-output <dir> \
        --broccoli-output <dir> --sonicparanoid2-output <dir> \
        --published-tsv orthobench_published.tsv \
        --replicates 1000 --sample-size 70 \
        --csv bootstrap_rank_support.csv

Requires benchmark.py to be importable (same directory by default, or set
--benchmark-path).
"""

import argparse
import csv
import importlib.util
import random
import sys
from collections import Counter, defaultdict
from pathlib import Path


TOOLS = ("orthofinder", "fastoma", "broccoli", "sonicparanoid2")

# Metric name -> (column header in orthobench_published.tsv, direction)
# direction "lower" means smaller values rank better; "higher" means larger
# values rank better.
METRICS = {
    "Missing_Refog": ("Missing RefOGs (%)", "lower"),
    "Missing_Genes": ("Missing Genes (%)", "lower"),
    "Fusion":        ("RefOG Fusions (%)", "lower"),
    "Fission":       ("RefOG Fissions (%)", "lower"),
    "Recall":        ("Recall", "higher"),
    "Precision":     ("Precision", "higher"),
    "Entropy":       ("Entropy", "lower"),
}

# Default row names in orthobench_published.tsv for each tool. Override on
# the command line if a different published table is used.
DEFAULT_ROWS = {
    "orthofinder": "OF3_Align_ST",
    "fastoma": "FastOMA",
    "broccoli": "Broccoli",
    "sonicparanoid2": "SP_def",
}


def load_benchmark_module(benchmark_path):
    spec = importlib.util.spec_from_file_location("benchmark", benchmark_path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def load_tool_predictions(bm, tool, tool_output):
    """Mirrors benchmark.py's main() tool-dispatch, minus the CLI/printing."""
    tool_output = Path(tool_output)
    if tool == "orthofinder":
        og_file = bm.find_orthofinder_og_file(tool_output)
        og_tbl = bm.load_orthofinder(og_file)
    elif tool == "fastoma":
        og_file = tool_output / "RootHOGs.tsv"
        if not og_file.exists():
            raise FileNotFoundError(f"FastOMA file not found: {og_file}")
        og_tbl = bm.load_fastoma(og_file)
    elif tool == "broccoli":
        og_file = tool_output / "dir_step3" / "orthologous_groups.txt"
        if not og_file.exists():
            raise FileNotFoundError(f"Broccoli file not found: {og_file}")
        og_tbl = bm.load_broccoli(og_file)
    elif tool == "sonicparanoid2":
        og_file = bm.find_sonicparanoid_og_file(tool_output)
        og_tbl = bm.load_sonicparanoid(og_file)
    else:
        raise ValueError(f"Unknown tool: {tool}")
    return bm.remove_singletons(og_tbl)


def index_by_refog(rows, refog_index):
    """rows: iterable of tuples whose element at refog_index is the refog id.
    Returns dict refog -> list of rows, for O(1) lookup per sampled refog."""
    index = defaultdict(list)
    for row in rows:
        index[row[refog_index]].append(row)
    return index


def load_published_values(published_tsv, row_map):
    values = {}
    with open(published_tsv, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        rows = {r["Methods"]: r for r in reader}
    for tool, row_name in row_map.items():
        if row_name not in rows:
            raise ValueError(
                f"Row '{row_name}' (for tool '{tool}') not found in {published_tsv}. "
                f"Available rows: {sorted(rows.keys())}"
            )
        values[tool] = rows[row_name]
    return values


def build_published_blocks(published_values, column, direction, tools):
    """Group tools into tie blocks based on published values, ordered best
    to worst. Returns a list of frozensets, e.g.
    [{"orthofinder","broccoli","sonicparanoid2"}, {"fastoma"}]."""
    vals = {t: float(published_values[t][column]) for t in tools}
    reverse = (direction == "higher")
    sorted_tools = sorted(tools, key=lambda t: vals[t], reverse=reverse)

    blocks = []
    current_block = [sorted_tools[0]]
    for t in sorted_tools[1:]:
        if abs(vals[t] - vals[current_block[-1]]) < 1e-9:
            current_block.append(t)
        else:
            blocks.append(frozenset(current_block))
            current_block = [t]
    blocks.append(frozenset(current_block))
    return blocks


def consistent_with_blocks(ordering, blocks):
    """ordering: tuple of tool names, best to worst (a full strict order for
    one replicate). blocks: published tie-block structure, best to worst.
    Returns True if `ordering` respects every between-block comparison in
    `blocks` (within-block order in `ordering` is unconstrained)."""
    position = {t: i for i, t in enumerate(ordering)}
    for i in range(len(blocks)):
        for j in range(i + 1, len(blocks)):
            for a in blocks[i]:
                for b in blocks[j]:
                    if position[a] >= position[b]:
                        return False
    return True


def rank_order(values_by_tool, direction, tools):
    """values_by_tool: dict tool -> metric value for one replicate.
    Returns a tuple of tool names, best to worst. Ties are broken by tool
    name for determinism (true numeric ties are effectively impossible with
    continuous simulated data, but this keeps the function well-defined)."""
    reverse = (direction == "higher")
    return tuple(sorted(tools, key=lambda t: (values_by_tool[t], t), reverse=reverse))


def published_pairwise_relations(blocks, tools):
    """blocks: published tie-block structure, best to worst. Returns dict
    {(a, b): relation} for every unordered pair a < b (alphabetically),
    where relation is "a_better", "b_better", or "tied" - derived purely
    from which block each tool falls in."""
    block_index = {t: i for i, block in enumerate(blocks) for t in block}
    relations = {}
    sorted_tools = sorted(tools)
    for i in range(len(sorted_tools)):
        for j in range(i + 1, len(sorted_tools)):
            a, b = sorted_tools[i], sorted_tools[j]
            if block_index[a] < block_index[b]:
                relations[(a, b)] = "a_better"
            elif block_index[a] > block_index[b]:
                relations[(a, b)] = "b_better"
            else:
                relations[(a, b)] = "tied"
    return relations


def replicate_pairwise_relations(ordering, pairs):
    """ordering: tuple of tool names, best to worst, for one replicate.
    pairs: iterable of (a, b) tuples (a < b alphabetically). Returns dict
    {(a, b): "a_better" or "b_better"} for this replicate (ties are
    effectively impossible with continuous simulated values)."""
    position = {t: i for i, t in enumerate(ordering)}
    return {
        (a, b): "a_better" if position[a] < position[b] else "b_better"
        for a, b in pairs
    }


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--output-folder", "-o", required=True,
                        help="Simulation output folder (contains Simulated_Orthogroups.txt).")
    parser.add_argument("--orthofinder-output", required=True)
    parser.add_argument("--fastoma-output", required=True)
    parser.add_argument("--broccoli-output", required=True)
    parser.add_argument("--sonicparanoid2-output", required=True)
    parser.add_argument("--published-tsv", required=True,
                        help="Path to orthobench_published.tsv.")
    parser.add_argument("--replicates", type=int, default=1000)
    parser.add_argument("--sample-size", type=int, default=70)
    parser.add_argument("--seed", type=int, default=None,
                        help="Optional seed for reproducibility. Omit for a "
                             "fresh nondeterministic draw each run (matches "
                             "the rest of the pipeline's convention).")
    parser.add_argument("--benchmark-path", default=None,
                        help="Path to benchmark.py. Defaults to a file named "
                             "benchmark.py in the same directory as this script.")
    parser.add_argument("--csv", default=None,
                        help="Optional path to write the results table as CSV.")
    parser.add_argument("--pairwise-csv", default=None,
                        help="Optional path to write the per-pair (6 pairs x 7 "
                             "metrics) rank-agreement breakdown as CSV.")
    parser.add_argument("--of3-row", default=DEFAULT_ROWS["orthofinder"])
    parser.add_argument("--fastoma-row", default=DEFAULT_ROWS["fastoma"])
    parser.add_argument("--broccoli-row", default=DEFAULT_ROWS["broccoli"])
    parser.add_argument("--sp2-row", default=DEFAULT_ROWS["sonicparanoid2"])

    args = parser.parse_args()

    benchmark_path = (
        Path(args.benchmark_path) if args.benchmark_path
        else Path(__file__).resolve().parent / "benchmark.py"
    )
    if not benchmark_path.exists():
        print(f"ERROR: could not find benchmark.py at {benchmark_path}. "
              "Pass --benchmark-path explicitly.")
        sys.exit(1)
    bm = load_benchmark_module(benchmark_path)

    row_map = {
        "orthofinder": args.of3_row,
        "fastoma": args.fastoma_row,
        "broccoli": args.broccoli_row,
        "sonicparanoid2": args.sp2_row,
    }

    output_folder = Path(args.output_folder).resolve()
    ref_file = output_folder / "Simulated_Orthogroups.txt"
    if not ref_file.exists():
        print(f"ERROR: reference file not found: {ref_file}")
        sys.exit(1)

    print(f"Simulation dir:  {output_folder}")
    print(f"Published table: {args.published_tsv}")
    print(f"Row mapping:     {row_map}")
    print(f"Replicates:      {args.replicates}")
    print(f"Sample size:     {args.sample_size}")
    print("")

    ref_og_ids, ref_genes = bm.load_reference(ref_file)
    ref_og_ids = set(ref_og_ids)
    all_refogs = sorted(ref_og_ids)
    if args.sample_size > len(all_refogs):
        print(f"ERROR: sample size {args.sample_size} exceeds available "
              f"reference orthogroups ({len(all_refogs)}).")
        sys.exit(1)

    ref_genes_index = index_by_refog(ref_genes, refog_index=0)

    tool_output_dirs = {
        "orthofinder": args.orthofinder_output,
        "fastoma": args.fastoma_output,
        "broccoli": args.broccoli_output,
        "sonicparanoid2": args.sonicparanoid2_output,
    }

    print("Loading tool predictions...")
    og_tbl_index = {}
    for tool in TOOLS:
        og_tbl = load_tool_predictions(bm, tool, tool_output_dirs[tool])
        og_tbl_index[tool] = index_by_refog(og_tbl, refog_index=2)
        print(f"  {tool}: {len(og_tbl)} gene rows loaded")
    print("")

    published_values = load_published_values(args.published_tsv, row_map)
    published_blocks = {
        metric: build_published_blocks(published_values, column, direction, TOOLS)
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

    if args.seed is not None:
        random.seed(args.seed)

    all_pairs = [(a, b) for i, a in enumerate(sorted(TOOLS))
                 for b in sorted(TOOLS)[i + 1:]]
    published_pairs = {
        metric: published_pairwise_relations(published_blocks[metric], TOOLS)
        for metric in METRICS
    }

    match_counts = {metric: 0 for metric in METRICS}
    ordering_counts = {metric: Counter() for metric in METRICS}
    pairwise_match_counts = {
        metric: {pair: 0 for pair in all_pairs} for metric in METRICS
    }

    for rep in range(args.replicates):
        sample = random.sample(all_refogs, args.sample_size)
        sample_set = set(sample)

        ref_genes_sample = [
            gene_row for refog in sample for gene_row in ref_genes_index.get(refog, [])
        ]

        og_tbl_sample = {}
        for tool in TOOLS:
            idx = og_tbl_index[tool]
            og_tbl_sample[tool] = [
                row for refog in sample for row in idx.get(refog, [])
            ]

        # Compute all seven metrics for each tool on this replicate's sample.
        values_by_metric = {metric: {} for metric in METRICS}
        for tool in TOOLS:
            missing_refog, missing_genes = bm.compute_missing(
                sample_set, ref_genes_sample, og_tbl_sample[tool]
            )
            fusion = bm.compute_fusion(sample_set, og_tbl_sample[tool])
            fission = bm.compute_fission(sample_set, og_tbl_sample[tool])
            recall, precision = bm.compute_recall_precision(
                sample_set, ref_genes_sample, og_tbl_sample[tool]
            )
            entropy = bm.compute_entropy(sample_set, ref_genes_sample, og_tbl_sample[tool])

            values_by_metric["Missing_Refog"][tool] = missing_refog
            values_by_metric["Missing_Genes"][tool] = missing_genes
            values_by_metric["Fusion"][tool] = fusion
            values_by_metric["Fission"][tool] = fission
            # benchmark.py's compute_recall_precision returns fractions
            # (0-1); orthobench_published.tsv reports these as percentages
            # (0-100), same as every other column - but since only rank
            # order is ever compared here (never magnitude), the scale
            # doesn't actually matter for correctness. Scaled anyway so
            # printed/exported values are self-explanatory if inspected.
            values_by_metric["Recall"][tool] = recall * 100
            values_by_metric["Precision"][tool] = precision * 100
            values_by_metric["Entropy"][tool] = entropy

        for metric, (_, direction) in METRICS.items():
            ordering = rank_order(values_by_metric[metric], direction, TOOLS)
            ordering_counts[metric][ordering] += 1
            if consistent_with_blocks(ordering, published_blocks[metric]):
                match_counts[metric] += 1

            rep_pairs = replicate_pairwise_relations(ordering, all_pairs)
            for pair in all_pairs:
                published_relation = published_pairs[metric][pair]
                if published_relation == "tied":
                    # Published data draws no distinction for this pair, so
                    # there is nothing to validate the resample against -
                    # counting every replicate as a "match" would be vacuous
                    # and misleadingly inflate support. Left untallied;
                    # reported as N/A below.
                    continue
                if rep_pairs[pair] == published_relation:
                    pairwise_match_counts[metric][pair] += 1

        if (rep + 1) % max(1, args.replicates // 10) == 0:
            print(f"  ...{rep + 1}/{args.replicates} replicates done")

    print("")
    print("RESULTS")
    print("")

    results_rows = []
    for metric in METRICS:
        n = args.replicates
        published_bucket_count = match_counts[metric]
        match_pct = published_bucket_count / n * 100

        # Build the set of entries to rank: the published-consistent bucket
        # as one entry, plus every distinct observed ordering that is NOT
        # consistent with the published order as its own entry.
        entries = [("<published order>", published_bucket_count)]
        for ordering, count in ordering_counts[metric].items():
            if not consistent_with_blocks(ordering, published_blocks[metric]):
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
        print(f"  Full published order reproduced in {published_bucket_count}/{n} "
              f"resamples ({match_pct:.1f}%)  [strict - all pairwise comparisons "
              f"at once; see pairwise breakdown below for where this comes from]")
        print(f"  Rank of published order by frequency: #{published_rank} "
              f"(out of {len(entries)} distinct patterns)")
        print(f"  Most frequent pattern overall: {top_str} "
              f"({top_count}/{n} = {top_count / n * 100:.1f}%)")
        print(f"  Distinct full orderings observed: {n_distinct_orderings}")
        print(f"  Pairwise breakdown:")
        for pair in all_pairs:
            a, b = pair
            relation = published_pairs[metric][pair]
            if relation == "tied":
                label = f"{a} = {b} (published tie)"
                print(f"    {label:35s} N/A - published data draws no "
                      f"distinction, nothing to validate")
            else:
                label = f"{a} > {b}" if relation == "a_better" else f"{b} > {a}"
                pair_count = pairwise_match_counts[metric][pair]
                pair_pct = pair_count / n * 100
                print(f"    {label:35s} reproduced in {pair_count}/{n} ({pair_pct:.1f}%)")
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

    if args.pairwise_csv:
        pairwise_rows = []
        for metric in METRICS:
            for pair in all_pairs:
                a, b = pair
                relation = published_pairs[metric][pair]
                if relation == "tied":
                    pairwise_rows.append({
                        "metric": metric,
                        "tool_a": a,
                        "tool_b": b,
                        "published_relation": relation,
                        "match_count": "NA",
                        "replicates": args.replicates,
                        "match_pct": "NA",
                    })
                else:
                    pair_count = pairwise_match_counts[metric][pair]
                    pairwise_rows.append({
                        "metric": metric,
                        "tool_a": a,
                        "tool_b": b,
                        "published_relation": relation,
                        "match_count": pair_count,
                        "replicates": args.replicates,
                        "match_pct": round(pair_count / args.replicates * 100, 2),
                    })
        pairwise_csv_path = Path(args.pairwise_csv)
        pairwise_csv_path.parent.mkdir(parents=True, exist_ok=True)
        with open(pairwise_csv_path, "w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=list(pairwise_rows[0].keys()))
            writer.writeheader()
            writer.writerows(pairwise_rows)
        print(f"Pairwise breakdown written to: {pairwise_csv_path}")


if __name__ == "__main__":
    main()
