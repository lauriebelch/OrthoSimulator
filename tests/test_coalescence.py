'''
Unit tests for the coalescent rearrangement used by the optional ILS step
(scripts/OrthogroupSimulation/coalescence.py).

The checks use small trees and known answers: closed-form coalescent
probabilities, the rule that copies from a duplication stay together, and
the direction of the theta effect (a large theta means little ILS).
'''

import math
import os
import random
import sys

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)

from ete4 import Tree

import scripts.OrthogroupSimulation.coalescence as co

TAG = "XTAGX"


def make_tree(newick):
    return Tree(newick, parser=1)


# a gene tree with a duplication (node 7) inside a larger tree.
# Names are <TAG><vertex id> and, like SaGePhy's, can have text after '_'.
# The tree is ultrametric (every tip is 3 from the root), as SaGePhy's
# pruned trees are. The rearrangement assumes this.
DUP_TREE = ("(((A1:1,A2:1)XTAGX7_7:1,B:2)XTAGX8_8:1,"
            "(C:2,D:2)XTAGX9_9:1)XTAGX10_10;")


def leaf_names(tree):
    return sorted(leaf.name for leaf in tree.leaves())


def sister_sets(tree):
    '''the set of leaf-name sets below every internal node'''
    return {frozenset(l.name for l in node.leaves())
            for node in tree.traverse() if not node.is_leaf}


##########################################
######### probabilities ##################
##########################################

def test_node_id_reads_the_vertex_number():
    assert co.node_id("XTAGX12_22", TAG) == "12"
    assert co.node_id("XTAGX3", TAG) == "3"


def test_coalesce_prob_one_lineage_is_certain():
    assert co.coalesce_all_prob(1, 0.0) == 1.0
    assert co.coalesce_all_prob(1, 5.0) == 1.0


def test_coalesce_prob_two_lineages_matches_closed_form():
    # two lineages coalesce at rate 1, so P(coalesced by u) = 1 - exp(-u)
    for u in (0.01, 0.5, 1.0, 3.0, 10.0):
        assert abs(co.coalesce_all_prob(2, u) - (1 - math.exp(-u))) < 1e-9


def test_coalesce_prob_rises_with_time_and_stays_between_0_and_1():
    for m in (2, 3, 5, 10):
        previous = -1.0
        for u in (0.05, 0.2, 1.0, 3.0, 8.0):
            p = co.coalesce_all_prob(m, u)
            assert 0.0 <= p <= 1.0
            assert p >= previous - 1e-12
            previous = p


def test_coalesce_prob_more_lineages_is_slower():
    u = 1.0
    assert co.coalesce_all_prob(2, u) > co.coalesce_all_prob(4, u) > co.coalesce_all_prob(8, u)


def test_float_and_high_precision_versions_agree_when_both_are_reliable():
    for m, u in ((2, 1.0), (3, 2.0), (5, 3.0), (8, 5.0)):
        assert abs(co.coalesce_all_prob(m, u) - co.coalesce_all_prob_mp(m, u)) < 1e-8


def test_many_lineages_in_a_short_window_does_not_break():
    # the float formula loses all precision here; the function must fall
    # back to high precision and still return a probability
    p = co.coalesce_all_prob(60, 0.05)
    assert 0.0 <= p <= 1.0 and math.isfinite(p)


def test_bounded_next_time_lies_inside_the_window():
    rng = random.Random(1)
    for k in (2, 3, 6):
        for T in (0.1, 1.0, 4.0):
            for _ in range(20):
                t = co.bounded_next_time(k, T, rng)
                assert 0.0 <= t <= T


##########################################
######### merging ########################
##########################################

def test_merge_two_joins_lineages_and_sets_branch_lengths():
    rng = random.Random(3)
    pool = [("a", 0.0), ("b", 0.0), ("c", 0.0)]
    co.merge_two(pool, 2.0, 4.0, rng)
    assert len(pool) == 2
    merged = [text for text, height in pool if text.startswith("(")][0]
    # two tips met at coalescent time 2.0, theta 4.0 -> each branch is 2.0 / 4.0
    assert merged.count(":0.5") == 2
    assert [height for text, height in pool if text == merged][0] == 2.0


def test_bounded_race_leaves_one_lineage():
    rng = random.Random(5)
    pool = [("a", 0.0), ("b", 0.0), ("c", 0.0), ("d", 0.0)]
    result = co.bounded_race(pool, 0.0, 3.0, 1.0, rng)
    assert len(result) == 1


##########################################
######### whole-tree behaviour ###########
##########################################

def test_rearranged_tree_keeps_the_same_leaves():
    tree = make_tree(DUP_TREE)
    for seed in range(10):
        new_tree, rf = co.main(tree, {"7"}, theta=2.0, seed=seed, tag=TAG)
        assert leaf_names(new_tree) == leaf_names(tree)
        assert 0.0 <= rf <= 1.0


def test_same_seed_gives_the_same_tree():
    tree = make_tree(DUP_TREE)
    a, _ = co.main(tree, {"7"}, theta=1.0, seed=42, tag=TAG)
    b, _ = co.main(tree, {"7"}, theta=1.0, seed=42, tag=TAG)
    assert a.write(parser=1) == b.write(parser=1)


def test_copies_from_a_duplication_always_stay_together():
    # A1 and A2 descend from the duplication at vertex 7, so they must be
    # sisters however strong the ILS is
    tree = make_tree(DUP_TREE)
    for seed in range(40):
        new_tree, _ = co.main(tree, {"7"}, theta=0.05, seed=seed, tag=TAG)
        assert frozenset({"A1", "A2"}) in sister_sets(new_tree), f"seed {seed}"


def test_rearranged_trees_are_binary():
    tree = make_tree(DUP_TREE)
    for seed in range(20):
        new_tree, _ = co.main(tree, {"7"}, theta=0.5, seed=seed, tag=TAG)
        for node in new_tree.traverse():
            assert node.is_leaf or len(node.children) == 2


def test_large_theta_means_almost_no_rearrangement():
    # a very large theta gives very long branches in coalescent units, so
    # lineages coalesce long before they can mix
    tree = make_tree(DUP_TREE)
    scores = [co.main(tree, {"7"}, theta=1e6, seed=s, tag=TAG)[1] for s in range(30)]
    assert sum(scores) / len(scores) < 0.02


def test_small_theta_rearranges_more_than_large_theta():
    tree = make_tree(DUP_TREE)
    small = [co.main(tree, {"7"}, theta=0.1, seed=s, tag=TAG)[1] for s in range(60)]
    large = [co.main(tree, {"7"}, theta=50.0, seed=s, tag=TAG)[1] for s in range(60)]
    assert sum(small) / len(small) > sum(large) / len(large)


def test_branch_lengths_are_positive():
    tree = make_tree(DUP_TREE)
    for seed in range(10):
        new_tree, _ = co.main(tree, {"7"}, theta=1.0, seed=seed, tag=TAG)
        for node in new_tree.traverse():
            if not node.is_root:
                assert node.dist > 0


def test_lineages_below_a_duplication_must_coalesce_before_it():
    # vertex 5 is a duplication whose two children are speciations (3 and 4).
    # Each side's two lineages must meet before the duplication, so
    # {s1a, s2a} and {s1b, s2b} are always clades however strong the ILS is
    newick = "((s1a:1,s2a:1)XTAGX3_3:1,(s1b:1,s2b:1)XTAGX4_4:1)XTAGX5_5;"
    tree = make_tree(newick)
    for seed in range(60):
        new_tree, _ = co.main(tree, {"5"}, theta=0.05, seed=seed, tag=TAG)
        clades = sister_sets(new_tree)
        assert frozenset({"s1a", "s2a"}) in clades, f"seed {seed}"
        assert frozenset({"s1b", "s2b"}) in clades, f"seed {seed}"


def test_without_a_duplication_the_same_lineages_can_mix_under_strong_ils():
    # the same tree with vertex 5 as a speciation lets the sides mix
    newick = "((s1a:1,s2a:1)XTAGX3_3:1,(s1b:1,s2b:1)XTAGX4_4:1)XTAGX5_5;"
    tree = make_tree(newick)
    mixed = 0
    for seed in range(60):
        new_tree, _ = co.main(tree, set(), theta=0.05, seed=seed, tag=TAG)
        clades = sister_sets(new_tree)
        if frozenset({"s1a", "s2a"}) not in clades or frozenset({"s1b", "s2b"}) not in clades:
            mixed += 1
    assert mixed > 0
