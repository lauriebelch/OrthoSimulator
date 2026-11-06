'''
Unit tests for how the true orthologs and duplications are read from the
SaGePhy files. These decide the ground truth, so they are checked on small
trees where every pair can be listed by hand.

A pair of genes is orthologous when their last common ancestor is a
speciation. If it is a duplication, the pair are paralogs.
'''

import os
import sys
import tempfile
from types import SimpleNamespace

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)

import scripts.OrthogroupSimulation.orthogroup_simulation_utils as og

G2H_HEADER = ("# GUEST-TO-HOST MAP\n"
              "Host tree:\t(sp1,sp2);\n"
              "Guest vertex name:\tGuest vertex ID:\tGuest vertex type:\tGuest vertex time:\n")


def make_run(tmp, newick, nodes, name="1"):
    '''Write <name>.pruned.tree and <name>.pruned.guest2host.
    nodes = [(vertex name, id, type)]'''
    temp = os.path.join(tmp, "temporary_files")
    os.makedirs(temp, exist_ok=True)
    with open(os.path.join(temp, f"{name}.pruned.tree"), "w") as f:
        f.write(newick)
    with open(os.path.join(temp, f"{name}.pruned.guest2host"), "w") as f:
        f.write(G2H_HEADER)
        for vertex, vid, kind in nodes:
            f.write(f"{vertex}\t{vid}\t{kind}\t0.5\t1\t1\n")
    return SimpleNamespace(o=tmp)


def read_pairs(tmp, name="1"):
    with open(os.path.join(tmp, "temporary_files", f"{name}.orthologpairs.txt")) as f:
        next(f)
        return {tuple(line.rstrip("\n").split("\t")) for line in f if line.strip()}


# a1 and a2 are copies from a duplication (vertex 3); b1 is the other species.
# The speciation (vertex 4) is the root.
TREE_ONE_DUP = "((a1:1,a2:1)XTAGX3_3:1,b1:2)XTAGX4_4;"
NODES_ONE_DUP = [("XTAGX3_3", 3, "DUPLICATION"), ("XTAGX4_4", 4, "SPECIATION")]


def test_pairs_across_a_speciation_are_orthologs_and_pairs_across_a_duplication_are_not():
    with tempfile.TemporaryDirectory() as tmp:
        args = make_run(tmp, TREE_ONE_DUP, NODES_ONE_DUP)
        og.Orthologs("1", args)
        pairs = read_pairs(tmp)
    assert pairs == {("a1", "b1"), ("a2", "b1")}
    assert ("a1", "a2") not in pairs          # paralogs


def test_a_tree_with_only_speciations_gives_every_cross_species_pair():
    newick = "((a:1,b:1)XTAGX3_3:1,(c:1,d:1)XTAGX4_4:1)XTAGX5_5;"
    nodes = [("XTAGX3_3", 3, "SPECIATION"), ("XTAGX4_4", 4, "SPECIATION"), ("XTAGX5_5", 5, "SPECIATION")]
    with tempfile.TemporaryDirectory() as tmp:
        args = make_run(tmp, newick, nodes)
        og.Orthologs("1", args)
        pairs = read_pairs(tmp)
    assert len(pairs) == 6                      # all pairs of 4 genes
    assert ("a", "b") in pairs and ("c", "d") in pairs


def test_orthologs_are_not_counted_twice_and_are_written_in_sorted_order():
    with tempfile.TemporaryDirectory() as tmp:
        args = make_run(tmp, TREE_ONE_DUP, NODES_ONE_DUP)
        og.Orthologs("1", args)
        with open(os.path.join(tmp, "temporary_files", "1.orthologpairs.txt")) as f:
            lines = f.read().strip().split("\n")[1:]
    assert len(lines) == len(set(lines))
    for line in lines:
        a, b = line.split("\t")
        assert a <= b


def test_a_tree_with_two_duplications_in_a_row():
    # ((a1,a2)D1,a3)D2 are all paralogs of each other; b1 is the ortholog of all three
    newick = "(((a1:1,a2:1)XTAGX3_3:1,a3:2)XTAGX4_4:1,b1:3)XTAGX5_5;"
    nodes = [("XTAGX3_3", 3, "DUPLICATION"), ("XTAGX4_4", 4, "DUPLICATION"), ("XTAGX5_5", 5, "SPECIATION")]
    with tempfile.TemporaryDirectory() as tmp:
        args = make_run(tmp, newick, nodes)
        og.Orthologs("1", args)
        pairs = read_pairs(tmp)
    assert pairs == {("a1", "b1"), ("a2", "b1"), ("a3", "b1")}


def test_a_polytomy_at_a_speciation_stops_the_run():
    newick = "(a:1,b:1,c:1)XTAGX3_3;"
    with tempfile.TemporaryDirectory() as tmp:
        args = make_run(tmp, newick, [("XTAGX3_3", 3, "SPECIATION")])
        try:
            og.Orthologs("1", args)
        except SystemExit:
            return
    raise AssertionError("a polytomy was accepted")


def test_sagephy_comments_in_the_tree_are_ignored():
    newick = ("((a1:1[ID=0 HOST=0],a2:1[ID=1])XTAGX3_3:1[ID=2 VERTEXTYPE=Duplication],"
              "b1:2[ID=3])XTAGX4_4:0[ID=4 DISCPT=(6,0)];")
    with tempfile.TemporaryDirectory() as tmp:
        args = make_run(tmp, newick, NODES_ONE_DUP)
        og.Orthologs("1", args)
        assert read_pairs(tmp) == {("a1", "b1"), ("a2", "b1")}


def test_ortholog_parsing_relabels_genes_and_adds_species():
    with tempfile.TemporaryDirectory() as tmp:
        temp = os.path.join(tmp, "temporary_files")
        os.makedirs(temp)
        with open(os.path.join(temp, "7.relabelmap.txt"), "w") as f:
            f.write("Old_Name\tNew_Name\nXTAGX0_0\tHomo_sapiens_7_1\nXTAGX1_1\tMus_musculus_7_1\n")
        with open(os.path.join(temp, "7.orthologpairs.txt"), "w") as f:
            f.write("Gene_1\tGene_2\nXTAGX0_0\tXTAGX1_1\n")
        og.OrthologParsing("7", SimpleNamespace(o=tmp))
        with open(os.path.join(temp, "7.orthologs.txt")) as f:
            rows = [line.rstrip("\n").split("\t") for line in f]
    assert rows[0] == ["Orthogroup", "Species_1", "Species_2", "Gene_1", "Gene_2"]
    assert rows[1] == ["7", "Homo_sapiens", "Mus_musculus", "Homo_sapiens_7_1", "Mus_musculus_7_1"]


def test_species_from_gene_handles_any_number_of_underscores_in_the_species_name():
    assert og.SpeciesFromGene("Ecoli_3_1") == "Ecoli"
    assert og.SpeciesFromGene("Homo_sapiens_3_12") == "Homo_sapiens"
    assert og.SpeciesFromGene("Candida_albicans_var_3_1") == "Candida_albicans_var"


def test_duplication_ids_are_read_from_the_guest_to_host_file():
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "x.guest2host")
        with open(path, "w") as f:
            f.write(G2H_HEADER)
            f.write("XTAGX3_3\t3\tDUPLICATION\t0.5\t1\t1\n"
                    "XTAGX4_4\t4\tSPECIATION\t0.9\t2\t2\n"
                    "XTAGX9_9\t9\tDUPLICATION\t0.2\t1\t1\n"
                    "XTAGX0_0\t0\tLEAF\t0.0\t0\t0\n")
        assert og.ReadDuplicationIDs(path) == {"3", "9"}


def test_a_tree_with_no_duplications_has_no_duplication_ids():
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "x.guest2host")
        with open(path, "w") as f:
            f.write(G2H_HEADER + "XTAGX4_4\t4\tSPECIATION\t0.9\t2\t2\n")
        assert og.ReadDuplicationIDs(path) == set()


def test_check_for_root_looks_for_a_discretisation_point_near_the_species_count():
    # a 10-species tree: an orthogroup counts as present at the root when some
    # node has DISCPT covering at least 90% of the species (9 or 10)
    species_tree = "(" + ",".join(f"s{i}:1" for i in range(10)) + ");"
    with tempfile.TemporaryDirectory() as tmp:
        os.makedirs(os.path.join(tmp, "temporary_files"))
        sp = os.path.join(tmp, "species.nwk")
        with open(sp, "w") as f:
            f.write(species_tree)
        args = SimpleNamespace(o=tmp, s=sp)
        for discpt, expected in (("(10,0)", 1), ("(9,0)", 1), ("(8,0)", 0), ("(3,1)", 0)):
            with open(os.path.join(tmp, "temporary_files", "1.unpruned.tree"), "w") as f:
                f.write(f"(a:1,b:1)XTAGX1_1:1[VERTEXTYPE=Speciation DISCPT={discpt}];")
            assert og.CheckForRoot("1", args) == expected, discpt
