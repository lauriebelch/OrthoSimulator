#!/usr/bin/env python3
'''
Helper for tests/run_example_checks.sh.

Each sub-command checks one thing about a finished run and prints PASS or FAIL
lines. The exit code is 1 if anything failed.

  sim                  check the output folder of an orthogroup simulation
  truth-to-orthofinder write the true orthogroups in OrthoFinder's file format
  perfect              check a benchmark CSV shows a perfect score
  theta-trend          check that more ILS (lower theta) gives more tree change
  compare              print how two simulation runs compare (information only)
'''

import argparse
import csv
import glob
import os
import re
import statistics
import sys

failures = 0


def ok(message):
    print("  PASS  " + message)


def bad(message):
    global failures
    failures += 1
    print("  FAIL  " + message)


def check(condition, message):
    if condition:
        ok(message)
    else:
        bad(message)
    return condition


##########################################
######### reading helpers ################
##########################################

def read_fasta(path):
    seqs = {}
    name = None
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                name = line[1:].strip()
                seqs[name] = []
            elif name is not None:
                seqs[name].append(line.strip())
    return {k: "".join(v) for k, v in seqs.items()}


def read_stats(folder):
    path = os.path.join(folder, "Simulated_Orthogroup_statistics.txt")
    with open(path) as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    return rows


def ils_scores(rows):
    # tree-change scores. 'NA' is recorded for trees too small to compare
    # (for example orthogroups with only two species), so it is skipped
    return [float(r["ils_norm_rf"]) for r in rows if r["ils_norm_rf"] not in ("", "NA")]


def orthogroup_of_gene(gene):
    # genes are named <species>_<orthogroup>_<copy>, and species can contain '_'
    return gene.rsplit("_", 2)[-2]


##########################################
######### sim ############################
##########################################

def cmd_sim(args):
    folder = args.dir
    print(f"Checking simulation output: {folder}")

    files_present = True
    for name in ("Simulated_Orthogroups.txt", "Simulated_orthologs.txt",
                 "Simulated_Orthogroup_statistics.txt", "Run_log.txt"):
        path = os.path.join(folder, name)
        if not check(os.path.isfile(path) and os.path.getsize(path) > 0, f"{name} exists and is not empty"):
            files_present = False
    if not files_present:
        print("  The run did not finish, so the remaining checks were skipped.")
        return

    alignments = sorted(glob.glob(os.path.join(folder, "alignment_files", "*.fa")))
    trees = sorted(glob.glob(os.path.join(folder, "tree_files", "*.tre")))
    check(len(alignments) >= 0.9 * args.n,
          f"alignments written for at least 90% of the {args.n} requested orthogroups (found {len(alignments)})")
    check(len(trees) == len(alignments),
          f"one tree per alignment ({len(trees)} trees, {len(alignments)} alignments)")

    # alignments: readable, rectangular, no all-gap sequence, gene names in the expected format
    problems = []
    for path in alignments:
        seqs = read_fasta(path)
        og = os.path.basename(path).split(".")[0]
        if len(seqs) < 2:
            problems.append(f"{os.path.basename(path)}: fewer than 2 sequences")
            continue
        if len({len(s) for s in seqs.values()}) != 1:
            problems.append(f"{os.path.basename(path)}: sequences differ in length")
        if any(set(s) <= {"-"} for s in seqs.values()):
            problems.append(f"{os.path.basename(path)}: a sequence is all gaps")
        for gene in seqs:
            if not re.match(r"^.+_\d+_\d+$", gene) or orthogroup_of_gene(gene) != og:
                problems.append(f"{os.path.basename(path)}: unexpected gene name {gene}")
                break
    check(not problems, "all alignments are readable, rectangular and gap-sequence free"
          + ("" if not problems else f" ({len(problems)} problems, first: {problems[0]})"))

    # true orthogroup table has one row per alignment
    with open(os.path.join(folder, "Simulated_Orthogroups.txt")) as f:
        header = f.readline().rstrip("\n").split("\t")
        rows = [line for line in f if line.strip()]
    check(header[0] == "Orthogroup" and len(header) > 2, "Simulated_Orthogroups.txt has the expected header")
    check(len(rows) == len(alignments),
          f"one row per alignment in Simulated_Orthogroups.txt ({len(rows)} rows)")

    # true ortholog pairs
    with open(os.path.join(folder, "Simulated_orthologs.txt")) as f:
        pair_header = f.readline().rstrip("\n").split("\t")
        pairs = [line.rstrip("\n").split("\t") for line in f if line.strip()]
    check(pair_header == ["Orthogroup", "Species_1", "Species_2", "Gene_1", "Gene_2"],
          "Simulated_orthologs.txt has the expected header")
    check(len(pairs) > 0, f"there are true ortholog pairs ({len(pairs)})")
    wrong = [p for p in pairs if len(p) != 5 or orthogroup_of_gene(p[3]) != p[0] or orthogroup_of_gene(p[4]) != p[0]]
    check(not wrong, "every ortholog pair lies inside one orthogroup")

    # statistics table and the ILS columns
    stats = read_stats(folder)
    check(len(stats) == len(alignments), f"one row per alignment in the statistics table ({len(stats)} rows)")
    check("ils_theta" in stats[0] and "ils_norm_rf" in stats[0], "statistics table has the ILS columns")

    temp = os.path.join(folder, "temporary_files")
    ils_trees = glob.glob(os.path.join(temp, "*.pruned_ils.tree"))
    original_trees = glob.glob(os.path.join(temp, "*.pruned.tree"))
    check(len(original_trees) >= len(alignments),
          "the original SaGePhy tree is kept for every orthogroup (the ground truth is read from it)")

    if args.ils == "off":
        check(all(r["ils_theta"] == "" and r["ils_norm_rf"] == "" for r in stats),
              "ILS off: no ILS values are recorded")
        check(len(ils_trees) == 0, "ILS off: no rearranged trees were made (the ILS step was skipped)")
    else:
        check(all(r["ils_theta"] != "" and abs(float(r["ils_theta"]) - args.theta) < 1e-9 for r in stats),
              f"ILS on: every orthogroup records theta = {args.theta}")
        values = ils_scores(stats)
        check(all(r["ils_norm_rf"] != "" for r in stats) and len(values) > 0 and all(0.0 <= v <= 1.0 for v in values),
              f"ILS on: tree-change scores are between 0 and 1 ({len(values)} of {len(stats)} orthogroups have one; "
              "the rest are too small to compare)")
        check(len(ils_trees) >= len(alignments), "ILS on: a rearranged tree was made for every orthogroup")
        if values:
            print(f"  INFO  mean tree-change score (normalised RF) = {statistics.mean(values):.3f}")


##########################################
######### truth-to-orthofinder ###########
##########################################

def cmd_truth_to_orthofinder(args):
    src = os.path.join(args.dir, "Simulated_Orthogroups.txt")
    out_dir = os.path.join(args.out, "Results_example", "Orthogroups")
    os.makedirs(out_dir, exist_ok=True)
    n = 0
    with open(src) as f, open(os.path.join(out_dir, "Orthogroups.txt"), "w") as w:
        f.readline()
        for line in f:
            parts = line.rstrip("\n").split("\t")
            genes = [g for cell in parts[1:] for g in cell.split(",") if g]
            w.write("OG%07d: %s\n" % (int(parts[0]), " ".join(genes)))
            n += 1
    print(f"wrote {n} orthogroups to {out_dir}/Orthogroups.txt")


##########################################
######### perfect ########################
##########################################

def cmd_perfect(args):
    with open(args.csv) as f:
        row = next(csv.DictReader(f))
    print(f"Checking benchmark scores in {args.csv}")
    for column in ("Recall", "Precision", "F1"):
        check(abs(float(row[column]) - 1.0) < 1e-6, f"{column} = {row[column]} (expected 1)")
    for column in ("Fusion", "Fission", "Missing_Genes", "Missing_Refog"):
        check(abs(float(row[column])) < 1e-6, f"{column} = {row[column]} (expected 0)")


##########################################
######### theta-trend ####################
##########################################

def cmd_theta_trend(args):
    means = {}
    for folder, theta in zip(args.dirs, args.thetas):
        values = ils_scores(read_stats(folder))
        check(len(values) > 0, f"theta = {theta}: tree-change scores were recorded ({len(values)})")
        if values:
            means[theta] = statistics.mean(values)
            print(f"  INFO  theta = {theta}: mean tree-change score = {means[theta]:.3f}")
    ordered = sorted(means)
    check(all(means[a] >= means[b] - 0.02 for a, b in zip(ordered, ordered[1:])),
          "lower theta gives more tree change (the score falls as theta rises)")
    if len(ordered) >= 2:
        check(means[ordered[0]] > means[ordered[-1]],
              f"the lowest theta ({ordered[0]}) changes trees more than the highest ({ordered[-1]})")


##########################################
######### compare ########################
##########################################

def cmd_compare(args):
    columns = ["num_genes", "num_species", "num_duplications", "num_losses", "seq_length",
               "total_branch_length"]
    print(f"A = {args.a}\nB = {args.b}")
    a_rows, b_rows = read_stats(args.a), read_stats(args.b)
    print(f"  INFO  orthogroups: A = {len(a_rows)}, B = {len(b_rows)}")
    for column in columns:
        def mean(rows):
            vals = [float(r[column]) for r in rows if r.get(column) not in (None, "")]
            return statistics.mean(vals) if vals else float("nan")
        ma, mb = mean(a_rows), mean(b_rows)
        rel = abs(ma - mb) / ma * 100 if ma else float("nan")
        print(f"  INFO  mean {column:<22} A = {ma:10.3f}   B = {mb:10.3f}   difference = {rel:5.1f}%")
    print("  (Runs are not seeded, so small differences are expected.)")


##########################################
######### main ###########################
##########################################

def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)

    p = sub.add_parser("sim")
    p.add_argument("--dir", required=True)
    p.add_argument("--n", type=int, required=True, help="number of orthogroups that were requested")
    p.add_argument("--ils", choices=["off", "on"], required=True)
    p.add_argument("--theta", type=float, default=None)
    p.set_defaults(func=cmd_sim)

    p = sub.add_parser("truth-to-orthofinder")
    p.add_argument("--dir", required=True)
    p.add_argument("--out", required=True)
    p.set_defaults(func=cmd_truth_to_orthofinder)

    p = sub.add_parser("perfect")
    p.add_argument("--csv", required=True)
    p.set_defaults(func=cmd_perfect)

    p = sub.add_parser("theta-trend")
    p.add_argument("--dirs", nargs="+", required=True)
    p.add_argument("--thetas", nargs="+", type=float, required=True)
    p.set_defaults(func=cmd_theta_trend)

    p = sub.add_parser("compare")
    p.add_argument("--a", required=True)
    p.add_argument("--b", required=True)
    p.set_defaults(func=cmd_compare)

    args = parser.parse_args()
    if args.command == "sim" and args.ils == "on" and args.theta is None:
        parser.error("--theta is required with --ils on")
    args.func(args)
    if failures:
        print(f"\n{failures} check(s) FAILED")
        sys.exit(1)


if __name__ == "__main__":
    main()
