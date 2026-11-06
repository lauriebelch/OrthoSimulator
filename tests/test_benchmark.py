'''
Unit tests for benchmarking/benchmark.py.

The metrics are checked on tiny orthogroup tables where the answer can be
worked out by hand: a perfect prediction, a split orthogroup (fission), a
merged pair (fusion), and missing genes. The file readers for each tool's
output format are checked on small files.
'''

import importlib.util
import math
import os
import sys
import tempfile
from pathlib import Path

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

spec = importlib.util.spec_from_file_location(
    "benchmark", os.path.join(REPO_ROOT, "benchmarking", "benchmark.py"))
bm = importlib.util.module_from_spec(spec)
spec.loader.exec_module(bm)


def write(path, text):
    with open(path, "w") as f:
        f.write(text)


# reference: two true orthogroups, 0 (four genes) and 1 (two genes)
REF_TABLE = ("Orthogroup\tSp_a\tSp_b\n"
             "0\tSp_a_0_1,Sp_a_0_2\tSp_b_0_1,Sp_b_0_2\n"
             "1\tSp_a_1_1\tSp_b_1_1\n")


def reference():
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "ref.txt")
        write(path, REF_TABLE)
        return bm.load_reference(path)


def table(*groups):
    '''groups = {orthogroup id: [gene, ...]} -> rows (og, gene, true og)'''
    rows = []
    for og_id, genes in groups:
        for g in genes:
            rows.append((og_id, g, int(g.split("_")[-2])))
    return rows


PERFECT = table(("OG1", ["Sp_a_0_1", "Sp_a_0_2", "Sp_b_0_1", "Sp_b_0_2"]),
                ("OG2", ["Sp_a_1_1", "Sp_b_1_1"]))


##########################################
######### reference ######################
##########################################

def test_reference_is_read_with_all_genes():
    ids, genes = reference()
    assert ids == {0, 1}
    assert sorted(genes) == sorted([(0, "Sp_a_0_1"), (0, "Sp_a_0_2"), (0, "Sp_b_0_1"),
                                    (0, "Sp_b_0_2"), (1, "Sp_a_1_1"), (1, "Sp_b_1_1")])


##########################################
######### metrics ########################
##########################################

def test_perfect_prediction_scores_perfectly():
    ids, genes = reference()
    assert bm.compute_missing(ids, genes, PERFECT) == (0.0, 0.0)
    assert bm.compute_fusion(ids, PERFECT) == 0.0
    assert bm.compute_fission(ids, PERFECT) == 0.0
    assert bm.compute_recall_precision(ids, genes, PERFECT) == (1.0, 1.0)
    assert bm.compute_entropy(ids, genes, PERFECT) == 0.0


def test_a_split_orthogroup_is_fission_and_lowers_recall_only():
    ids, genes = reference()
    split = table(("OG1", ["Sp_a_0_1", "Sp_a_0_2"]), ("OG3", ["Sp_b_0_1", "Sp_b_0_2"]),
                  ("OG2", ["Sp_a_1_1", "Sp_b_1_1"]))
    assert bm.compute_fission(ids, split) == 50.0          # one of two true OGs is split
    assert bm.compute_fusion(ids, split) == 0.0
    recall, precision = bm.compute_recall_precision(ids, genes, split)
    assert precision == 1.0
    # true OG 0 (4 genes) has each predicted piece covering half of it: recall 0.5
    # true OG 1 is intact: recall 1.0. Weighted by size: (4*0.5 + 2*1.0) / 6
    assert abs(recall - (4 * 0.5 + 2 * 1.0) / 6) < 1e-12


def test_merging_two_true_orthogroups_is_fusion_and_lowers_precision_only():
    ids, genes = reference()
    fused = table(("OG1", ["Sp_a_0_1", "Sp_a_0_2", "Sp_b_0_1", "Sp_b_0_2", "Sp_a_1_1", "Sp_b_1_1"]))
    assert bm.compute_fusion(ids, fused) == 100.0          # both true OGs are in a fused group
    assert bm.compute_fission(ids, fused) == 0.0
    recall, precision = bm.compute_recall_precision(ids, genes, fused)
    assert recall == 1.0
    # predicted group holds 6 genes; 4 belong to OG 0 and 2 to OG 1
    assert abs(precision - (4 * (4 / 6) + 2 * (2 / 6)) / 6) < 1e-12


def test_missing_genes_and_missing_orthogroups_are_counted():
    ids, genes = reference()
    partial = table(("OG1", ["Sp_a_0_1", "Sp_a_0_2", "Sp_b_0_1"]))   # OG 1 absent, one gene of OG 0 absent
    missing_refog, missing_genes = bm.compute_missing(ids, genes, partial)
    assert missing_refog == 50.0
    assert abs(missing_genes - 3 / 6 * 100) < 1e-12


def test_entropy_is_zero_when_intact_and_positive_when_split():
    ids, genes = reference()
    split = table(("OG1", ["Sp_a_0_1", "Sp_a_0_2"]), ("OG3", ["Sp_b_0_1", "Sp_b_0_2"]),
                  ("OG2", ["Sp_a_1_1", "Sp_b_1_1"]))
    # true OG 0 (4 genes) is split into two halves: entropy 1 bit, divided by
    # log2(4) = 2 gives 0.5. True OG 1 contributes 0. Weighted by size: 4 * 0.5 / 6
    assert abs(bm.compute_entropy(ids, genes, split) - 4 * 0.5 / 6) < 1e-12


def test_all_metrics_stay_in_range_for_a_scrambled_prediction():
    ids, genes = reference()
    scrambled = table(("OG1", ["Sp_a_0_1", "Sp_b_1_1"]), ("OG2", ["Sp_a_0_2", "Sp_a_1_1", "Sp_b_0_1"]),
                      ("OG3", ["Sp_b_0_2"]))
    scrambled = bm.remove_singletons(scrambled)
    recall, precision = bm.compute_recall_precision(ids, genes, scrambled)
    assert 0.0 <= recall <= 1.0 and 0.0 <= precision <= 1.0
    assert 0.0 <= bm.compute_entropy(ids, genes, scrambled) <= 1.0
    assert 0.0 <= bm.compute_fusion(ids, scrambled) <= 100.0
    assert 0.0 <= bm.compute_fission(ids, scrambled) <= 100.0


def test_singleton_predicted_groups_are_dropped():
    rows = table(("OG1", ["Sp_a_0_1"]), ("OG2", ["Sp_a_1_1", "Sp_b_1_1"]))
    assert {r[0] for r in bm.remove_singletons(rows)} == {"OG2"}


##########################################
######### file readers ###################
##########################################

def test_orthofinder_reader():
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "Orthogroups.txt")
        write(path, "OG0000000: Sp_a_0_1 Sp_a_0_2 Sp_b_0_1\nOG0000001: Sp_a_1_1 Sp_b_1_1\n\n")
        rows = bm.load_orthofinder(path)
    assert ("OG0000000", "Sp_a_0_2", 0) in rows
    assert ("OG0000001", "Sp_b_1_1", 1) in rows
    assert len(rows) == 5


def test_orthofinder_reader_drops_protein_suffix_and_skips_odd_names():
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "Orthogroups.txt")
        write(path, "OG0: Sp_a_0_1.1 weird Sp_b_x_1\n")
        rows = bm.load_orthofinder(path)
    assert rows == [("OG0", "Sp_a_0_1", 0)]


def test_fastoma_reader_skips_the_header_row():
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "RootHOGs.tsv")
        write(path, "RootHOG\tProtein\nHOG:0001\tSp_a_0_1\nHOG:0001\tSp_b_0_1\nHOG:0002\tSp_a_1_1\n")
        rows = bm.load_fastoma(path)
    assert [r[1] for r in rows] == ["Sp_a_0_1", "Sp_b_0_1", "Sp_a_1_1"]
    assert rows[0][0] == "HOG:0001"


def test_broccoli_reader_skips_header_and_stars():
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "orthologous_groups.txt")
        write(path, "#OG_name\tProtein_names\nOG_1\tSp_a_0_1 Sp_b_0_1 *\nOG_2\tSp_a_1_1\n")
        rows = bm.load_broccoli(path)
    assert [(r[0], r[1]) for r in rows] == [("OG_1", "Sp_a_0_1"), ("OG_1", "Sp_b_0_1"), ("OG_2", "Sp_a_1_1")]


def test_sonicparanoid_reader_reads_every_species_column():
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "flat.ortholog_groups.tsv")
        write(path, "group_id\tSp_a\tSp_b\n1\tSp_a_0_1,Sp_a_0_2\tSp_b_0_1\n2\t*\tSp_b_1_1\n")
        rows = bm.load_sonicparanoid(path)
    assert sorted(r[1] for r in rows) == ["Sp_a_0_1", "Sp_a_0_2", "Sp_b_0_1", "Sp_b_1_1"]


def test_sonicparanoid_reader_needs_a_group_id_column():
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "x.tsv")
        write(path, "id\tSp_a\n1\tSp_a_0_1\n")
        try:
            bm.load_sonicparanoid(path)
        except ValueError:
            return
    raise AssertionError("a table without group_id was accepted")


def test_species_names_with_any_number_of_underscores_are_read_correctly():
    # simulated genes are <species>_<orthogroup>_<copy>, and species names
    # can have any number of underscores
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "Orthogroups.txt")
        write(path, "OG0: Ecoli_7_1 Homo_sapiens_7_2 Candida_albicans_var_7_1\n")
        rows = bm.load_orthofinder(path)
    assert [r[2] for r in rows] == [7, 7, 7]


def test_orthofinder_file_is_found_in_a_results_folder():
    with tempfile.TemporaryDirectory() as tmp:
        folder = Path(tmp) / "Results_Nov06" / "Orthogroups"
        folder.mkdir(parents=True)
        write(folder / "Orthogroups.txt", "OG0: a_0_1\n")
        assert bm.find_orthofinder_og_file(Path(tmp)) == folder / "Orthogroups.txt"


def test_missing_orthofinder_file_is_a_clear_error():
    with tempfile.TemporaryDirectory() as tmp:
        try:
            bm.find_orthofinder_og_file(Path(tmp))
        except FileNotFoundError as err:
            assert "Orthogroups.txt" in str(err)
            return
    raise AssertionError("no error for a missing file")
