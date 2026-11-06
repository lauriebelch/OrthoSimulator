'''
Unit tests for the helpers used by --Get-Parameters:
the trimmed summary statistics, the gap-position profile of an alignment,
the species and gene counts of an alignment, and species-overlap
duplication counting on gene trees.
'''

import math
import os
import sys
import tempfile

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)

from ete4 import Tree

import scripts.GetParameters.EstimateParameters as ep
import scripts.GetParameters.CallDuplications as cd


def write(path, text):
    with open(path, "w") as f:
        f.write(text)


##########################################
######### trimmed summaries ##############
##########################################

def test_trim_and_summarize_of_nothing_returns_nones():
    assert ep.trim_and_summarize([]) == (None, None, None, None)


def test_trimming_removes_extreme_values_from_the_mean():
    data = [1.0] * 8 + [100.0, -100.0]     # one extreme value each side
    mean, sd, _, _ = ep.trim_and_summarize(data)
    assert abs(mean - 1.0) < 1e-9
    assert sd == 0.0


def test_robust_sd_uses_the_median_absolute_deviation():
    data = list(range(1, 101))
    mean, sd, log_mean, log_sd = ep.trim_and_summarize(data)
    assert abs(mean - 50.5) < 1.0
    # the 10-90 percentile trim keeps about 10.9 to 90.1, so MAD is about 20
    assert 25 < sd < 33
    assert abs(log_mean - math.log(mean)) < 1e-12
    assert abs(log_sd - math.sqrt(math.log(1 + (sd / mean) ** 2))) < 1e-12


def test_log_values_are_none_for_a_nonpositive_mean():
    mean, sd, log_mean, log_sd = ep.trim_and_summarize([-3.0] * 10)
    assert log_mean is None and log_sd is None


def test_safe_format():
    assert ep.safe_format(None) == "N/A"
    assert ep.safe_format(1.23456, 2) == "1.23"


##########################################
######### alignment helpers ##############
##########################################

def test_species_and_gene_counts_from_headers():
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "og.fa")
        write(path, ">1_geneA\nMKV\n>1_geneB\nMKV\n>2_geneC\nMKV\n")
        assert ep.get_alignment_species_and_gene_counts(path) == (2, 3)


def test_fasta_sequences_are_counted():
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "og.fa")
        write(path, ">a\nMK\n>b\nMK\n>c\nMK\n")
        assert ep.count_fasta_seqs(path) == 3


def test_gap_profile_puts_gaps_in_the_right_percent_bin():
    # 100 columns, gaps (in 2 sequences) at column 50 -> bin 50
    seq1 = "A" * 50 + "-" + "A" * 49
    seq2 = "A" * 50 + "-" + "A" * 49
    seq3 = "A" * 100
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "a.fa")
        write(path, f">1\n{seq1}\n>2\n{seq2}\n>3\n{seq3}\n")
        counts, probs = ep.get_alignment_gap_profile(path)
    assert counts[50] == 2
    assert sum(counts.values()) == 2
    assert probs[50] == 1.0
    assert abs(sum(probs.values()) - 1.0) < 1e-12


def test_gap_profile_of_a_gapless_alignment_is_all_zeros():
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "a.fa")
        write(path, ">1\nAAAA\n>2\nAAAA\n")
        counts, probs = ep.get_alignment_gap_profile(path)
    assert sum(counts.values()) == 0 and sum(probs.values()) == 0.0


def test_gap_profile_has_exactly_100_bins_and_the_last_column_goes_in_bin_99():
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "a.fa")
        write(path, ">1\n" + "A" * 9 + "-\n>2\n" + "A" * 10 + "\n")
        counts, _ = ep.get_alignment_gap_profile(path)
    assert len(counts) == 100
    assert counts[90] == 1          # column 9 of 10 sits at 90%


def test_alignment_reader_matches_the_simulator_reader():
    import scripts.OrthogroupSimulation.orthogroup_simulation_utils as og
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "a.fa")
        write(path, ">x\nMK-\nLLA\n>y\nMKVLLA\n")
        assert ep.read_fasta_alignment(path) == og.ReadAlignment(path)


##########################################
######### duplications ###################
##########################################

def test_species_id_is_the_part_before_the_first_underscore():
    assert cd.GetSpeciesID("12_geneB") == "12"
    assert cd.GetSpeciesID("3_a_b_c") == "3"


def tree(newick):
    return Tree(newick, parser=1)


def test_a_tree_with_one_gene_per_species_has_no_duplications():
    assert cd.CountSpeciesOverlapDuplications(tree("((1_a,2_a),(3_a,4_a));")) == 0


def test_a_species_appearing_on_both_sides_of_a_node_is_a_duplication():
    assert cd.CountSpeciesOverlapDuplications(tree("((1_a,1_b),(2_a,3_a));")) == 1


def test_nested_duplications_are_each_counted():
    # (1_a,1_b) is one duplication; joining it with 1_c is a second
    assert cd.CountSpeciesOverlapDuplications(tree("(((1_a,1_b),1_c),2_a);")) == 2


def test_species_in_a_tree_are_listed_once():
    t = tree("((1_a,1_b),(2_a,3_a));")
    assert cd.GetRepresentedSpecies(t) == {"1", "2", "3"}


def test_genes_are_picked_out_by_species():
    t = tree("((1_a,1_b),(2_a,3_a));")
    assert sorted(cd.GetGenesFromSpecies(t, {"1", "3"})) == ["1_a", "1_b", "3_a"]
