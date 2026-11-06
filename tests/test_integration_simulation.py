'''
Integration tests: run the real simulator on the example data (a few
orthogroups, 2 threads) and check that the files it writes agree with each
other.

The runs are made once and shared by the tests. They need Java and the
bundled IQ-TREE and SaGePhy. If those cannot run on this machine the tests
print a note and pass without checking anything.
'''

import atexit
import glob
import os
import shutil
import subprocess
import sys
import tempfile
from types import SimpleNamespace

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TESTS_DIR = os.path.join(REPO_ROOT, "tests")
for p in (REPO_ROOT, TESTS_DIR):
    if p not in sys.path:
        sys.path.insert(0, p)

import check_outputs as co
import scripts.OrthogroupSimulation.orthogroup_simulation_utils as og

N_OG = 6
CONFIG = os.path.join("ExampleData", "example_simulation_config.txt")


def can_run():
    if shutil.which("java") is None:
        return False
    iqtree = os.path.join(REPO_ROOT, "bin", "iqtree3")
    return os.access(iqtree, os.X_OK)


_RUNS = {}
_WORKDIR = tempfile.mkdtemp(prefix="orthosim_integration_")
atexit.register(shutil.rmtree, _WORKDIR, ignore_errors=True)


def run_simulator(name, *extra, n=N_OG):
    '''run once per name; returns (output folder, finished process)'''
    if name not in _RUNS:
        out = os.path.join(_WORKDIR, name)
        proc = subprocess.run(
            [sys.executable, "-B", "OrthoSimulator.py", "--Orthogroup-simulation",
             "--output", out, "--Config", CONFIG, "--Orthogroups", str(n), "--threads", "2", *extra],
            cwd=REPO_ROOT, capture_output=True, text=True, timeout=900)
        _RUNS[name] = (out, proc)
    return _RUNS[name]


def skip_note(name):
    print(f"      (skipped {name}: Java or the bundled IQ-TREE cannot run here)")


def sim_check(folder, n, ils, theta=None):
    co.failures = 0
    co.cmd_sim(SimpleNamespace(dir=folder, n=n, ils=ils, theta=theta))
    assert co.failures == 0, f"{co.failures} output checks failed in {folder}"


def alignment_names(folder, og_id):
    names, seqs = og.ReadAlignment(os.path.join(folder, "alignment_files", f"{og_id}.alignment.fa"))
    return names, seqs


##########################################
######### run without ILS ################
##########################################

def test_run_without_ils_passes_the_output_checks():
    if not can_run():
        return skip_note("test_run_without_ils_passes_the_output_checks")
    out, proc = run_simulator("off")
    assert proc.returncode == 0, proc.stdout[-2000:] + proc.stderr[-2000:]
    sim_check(out, N_OG, "off")


def test_trees_and_alignments_have_the_same_genes():
    if not can_run():
        return skip_note("test_trees_and_alignments_have_the_same_genes")
    out, _ = run_simulator("off")
    from ete4 import Tree
    for tree_path in glob.glob(os.path.join(out, "tree_files", "*.tre")):
        og_id = os.path.basename(tree_path).split(".")[0]
        names, _ = alignment_names(out, og_id)
        tree = Tree(open(tree_path).read().strip(), parser=1)
        assert sorted(l.name for l in tree.leaves()) == sorted(names), og_id


def test_orthogroup_table_lists_exactly_the_genes_in_the_alignments():
    if not can_run():
        return skip_note("test_orthogroup_table_lists_exactly_the_genes_in_the_alignments")
    out, _ = run_simulator("off")
    with open(os.path.join(out, "Simulated_Orthogroups.txt")) as f:
        header = f.readline().rstrip("\n").split("\t")
        for line in f:
            parts = line.rstrip("\n").split("\t")
            genes = sorted(g for cell in parts[1:] for g in cell.split(",") if g)
            names, _ = alignment_names(out, parts[0])
            assert genes == sorted(names)
            # each gene sits in the column of its own species
            for col, cell in zip(header[1:], parts[1:]):
                for g in filter(None, cell.split(",")):
                    assert og.SpeciesFromGene(g) == col


def test_true_orthologs_can_be_recomputed_from_the_saved_sagephy_trees():
    if not can_run():
        return skip_note("test_true_orthologs_can_be_recomputed_from_the_saved_sagephy_trees")
    out, _ = run_simulator("off")
    copy = os.path.join(_WORKDIR, "off_copy")
    if not os.path.exists(copy):
        shutil.copytree(out, copy)
    with open(os.path.join(out, "Simulated_orthologs.txt")) as f:
        next(f)
        reported = {}
        for line in f:
            og_id, s1, s2, g1, g2 = line.rstrip("\n").split("\t")
            reported.setdefault(og_id, set()).add(tuple(sorted((g1, g2))))
    for tree_path in glob.glob(os.path.join(copy, "temporary_files", "*.pruned.tree")):
        og_id = os.path.basename(tree_path).split(".")[0]
        if not os.path.exists(os.path.join(out, "alignment_files", f"{og_id}.alignment.fa")):
            continue
        args = SimpleNamespace(o=copy)
        og.Orthologs(og_id, args)
        og.OrthologParsing(og_id, args)
        with open(os.path.join(copy, "temporary_files", f"{og_id}.orthologs.txt")) as f:
            next(f)
            again = {tuple(sorted(line.rstrip("\n").split("\t")[3:5])) for line in f if line.strip()}
        assert again == reported.get(og_id, set()), og_id


def test_ortholog_pairs_are_always_between_different_species():
    if not can_run():
        return skip_note("test_ortholog_pairs_are_always_between_different_species")
    out, _ = run_simulator("off")
    with open(os.path.join(out, "Simulated_orthologs.txt")) as f:
        next(f)
        for line in f:
            _, s1, s2, g1, g2 = line.rstrip("\n").split("\t")
            assert s1 != s2 and og.SpeciesFromGene(g1) == s1 and og.SpeciesFromGene(g2) == s2


def test_statistics_agree_with_the_alignments():
    if not can_run():
        return skip_note("test_statistics_agree_with_the_alignments")
    out, _ = run_simulator("off")
    for row in co.read_stats(out):
        names, seqs = alignment_names(out, row["orthogroup"])
        assert int(row["num_genes"]) == len(names)
        assert int(row["num_species"]) == len({og.SpeciesFromGene(n) for n in names})


def test_a_finished_output_folder_is_not_overwritten():
    if not can_run():
        return skip_note("test_a_finished_output_folder_is_not_overwritten")
    out, _ = run_simulator("off")
    before = sorted(os.listdir(out))
    proc = subprocess.run(
        [sys.executable, "-B", "OrthoSimulator.py", "--Orthogroup-simulation", "--output", out,
         "--Config", CONFIG, "--Orthogroups", "2", "--threads", "1"],
        cwd=REPO_ROOT, capture_output=True, text=True, timeout=300)
    assert "already exists" in proc.stdout
    assert sorted(os.listdir(out)) == before


def test_a_species_without_pfam_data_gives_a_clear_error():
    if not can_run():
        return skip_note("test_a_species_without_pfam_data_gives_a_clear_error")
    out = os.path.join(_WORKDIR, "no_pfam")
    proc = subprocess.run(
        [sys.executable, "-B", "OrthoSimulator.py", "--Orthogroup-simulation", "--output", out,
         "--Config", CONFIG, "--Orthogroups", "2", "--threads", "1", "--species", "Not_a_species"],
        cwd=REPO_ROOT, capture_output=True, text=True, timeout=300)
    assert proc.returncode != 0
    assert "PFAM" in proc.stderr and "Not_a_species" in proc.stderr


##########################################
######### run with ILS ###################
##########################################

def test_run_with_ils_passes_the_output_checks():
    if not can_run():
        return skip_note("test_run_with_ils_passes_the_output_checks")
    out, proc = run_simulator("on", "--ils-theta", "3")
    assert proc.returncode == 0, proc.stdout[-2000:] + proc.stderr[-2000:]
    sim_check(out, N_OG, "on", theta=3.0)


def test_ils_does_not_change_which_genes_are_in_each_tree():
    if not can_run():
        return skip_note("test_ils_does_not_change_which_genes_are_in_each_tree")
    out, _ = run_simulator("on", "--ils-theta", "3")
    from ete4 import Tree
    for tree_path in glob.glob(os.path.join(out, "tree_files", "*.tre")):
        og_id = os.path.basename(tree_path).split(".")[0]
        names, _ = alignment_names(out, og_id)
        tree = Tree(open(tree_path).read().strip(), parser=1)
        assert sorted(l.name for l in tree.leaves()) == sorted(names)


def test_with_ils_the_truth_still_comes_from_the_original_tree():
    if not can_run():
        return skip_note("test_with_ils_the_truth_still_comes_from_the_original_tree")
    out, _ = run_simulator("on", "--ils-theta", "3")
    copy = os.path.join(_WORKDIR, "on_copy")
    if not os.path.exists(copy):
        shutil.copytree(out, copy)
    with open(os.path.join(out, "Simulated_orthologs.txt")) as f:
        next(f)
        reported = {}
        for line in f:
            og_id, s1, s2, g1, g2 = line.rstrip("\n").split("\t")
            reported.setdefault(og_id, set()).add(tuple(sorted((g1, g2))))
    checked = 0
    for tree_path in glob.glob(os.path.join(copy, "temporary_files", "*.pruned.tree")):
        og_id = os.path.basename(tree_path).split(".")[0]
        if not os.path.exists(os.path.join(out, "alignment_files", f"{og_id}.alignment.fa")):
            continue
        args = SimpleNamespace(o=copy)
        og.Orthologs(og_id, args)
        og.OrthologParsing(og_id, args)
        with open(os.path.join(copy, "temporary_files", f"{og_id}.orthologs.txt")) as f:
            next(f)
            again = {tuple(sorted(line.rstrip("\n").split("\t")[3:5])) for line in f if line.strip()}
        assert again == reported.get(og_id, set())
        checked += 1
    assert checked > 0


def test_zero_theta_is_refused_before_any_simulation_starts():
    if not can_run():
        return skip_note("test_zero_theta_is_refused_before_any_simulation_starts")
    out = os.path.join(_WORKDIR, "theta0")
    proc = subprocess.run(
        [sys.executable, "-B", "OrthoSimulator.py", "--Orthogroup-simulation", "--output", out,
         "--Config", CONFIG, "--Orthogroups", "2", "--threads", "1", "--ils-theta", "0"],
        cwd=REPO_ROOT, capture_output=True, text=True, timeout=300)
    assert proc.returncode != 0
    assert "ils_theta must be > 0" in proc.stderr
    assert not glob.glob(os.path.join(out, "alignment_files", "*"))
