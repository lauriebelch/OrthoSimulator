#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Post-hoc multispecies-coalescent (MSC) rearrangement of a simulated gene
tree, introducing incomplete lineage sorting (ILS) at a severity controlled
by a single parameter, theta (tree depth in coalescent units: larger =
less ILS, smaller = more ILS).

The OrthoSimulator paper (ILS methods) gives the model and references. In
short: SaGePhy's GuestTreeGen assumes gene lineages coalesce instantaneously
at each speciation, so it cannot produce ILS on its own. This module
rearranges the simulated topology to be consistent with the MSC, treating
duplication events as fixed, dated ground truth (via the bounded coalescent),
while leaving losses untouched (the pruned tree already has no signal from
lost lineages, so ignoring them is exact, not an approximation).

This module is intentionally free of any pipeline-specific I/O (no file
paths, no SaGePhy subprocess calls, no logging) — the orthogroup simulation
pipeline (orthogroup_simulation_utils.py) is responsible for reading the
input tree/duplication IDs and writing the output tree.
"""

import random
import math
from ete4 import Tree
from mpmath import mp, mpf, exp as mpexp

################################################
## core functions
################################################

## function that will rearrange a gene tree topology by doing coalescence
## backwards in time and defining a new tree based on where lineages actually
## coalesce
def rearrange(tree, theta, seed, dup_ids, tag):
    rng = random.Random(seed)
    # measure tree height down one path
    # theta / height turns a branch length into an amount of coalescent time
    height = 0.0
    node = tree
    while len(node.children) > 0:
        node = node.children[0]
        height = height + node.dist
    coalescent_units = theta / height
    # ultimately we track lineages 'entering' the branch from the tip end
    # and then see if they coalesce in this branch, or get passed 'up'
    # where they might cause ILS

    # for each node remember the lineages leaving the top of its branch, and the
    # coalescent time there. A lineage is a pair (some newick_text, birth_height).
    lineages_at_top = {}
    time_at_top = {}
    for node in tree.traverse("postorder"):
        # gather the lineages entering the bottom of this node's branch
        if len(node.children) == 0:
            pool = [(node.name, 0.0)]
            start_time = 0.0
        else:
            pool = []
            for child in node.children:
                pool = pool + lineages_at_top[child]
            start_time = time_at_top[node.children[0]]

        # Duplications are annoying because we have to use a bounded coalescent
        # to ensure that the duplication lineages have coalesced by the time we
        # get to the duplication node (and dont ALL coalesce there, causing
        # a billion unwanted polytomies)

        # duplication children arrive as single lineages.
        # here they meet and merge at the duplication point itself
        if len(node.children) > 0 and node_id(node.name, tag) in dup_ids:
            while len(pool) > 1:
                merge_two(pool, start_time, theta, rng)
        # if we are at the root, then there are no branches above where lineages
        # could go. Therefore, we force whatever is left to merge
        if node.up is None:
            time = start_time
            while len(pool) > 1:
                k = len(pool)
                time = time + rng.expovariate(k * (k - 1) / 2.0)
                merge_two(pool, time, theta, rng)
            return pool[0][0] + ";"
        top_time = start_time + coalescent_units * node.dist
        # a branch feeding a duplication must fully coalesce BEFORE the boundary
        # because the duplication event is the truth of where the duplication happened
        if node_id(node.up.name, tag) in dup_ids:
            pool = bounded_race(pool, start_time, top_time, theta, rng)
        # an ordinary branch = race, and carry any leftovers up
        else:
            time = start_time
            while len(pool) > 1 and time < top_time:
                k = len(pool)
                wait = rng.expovariate(k * (k - 1) / 2.0)
                if time + wait > top_time:
                    break
                time = time + wait
                merge_two(pool, time, theta, rng)
        lineages_at_top[node] = pool
        time_at_top[node] = top_time


## function that picks two lineages to merge
def merge_two(pool, time, theta, rng):
    a_text, a_height = pool.pop(rng.randrange(len(pool)))
    b_text, b_height = pool.pop(rng.randrange(len(pool)))
    a_length = (time - a_height) / theta
    b_length = (time - b_height) / theta
    combined = "(" + a_text + ":" + str(round(a_length, 10)) + \
               "," + b_text + ":" + str(round(b_length, 10)) + ")"
    pool.append((combined, time))


## probability of coalescence
## fast float version, but it suffers catastrophic cancellation when there are
## many lineages coalescing in a short window (short branch / hundreds of
## species), so in that case we fall back to arbitrary precision (mpmath)
def coalesce_all_prob(m, u):
    if m <= 1:
        return 1.0
    rates = [j * (j - 1) / 2.0 for j in range(2, m + 1)]
    total = 0.0
    maxterm = 0.0
    overflow = False
    for j in range(len(rates)):
        term = 1.0
        for i in range(len(rates)):
            if i != j:
                term = term * rates[i] / (rates[i] - rates[j])
        term = term * math.exp(-rates[j] * u)
        total = total + term
        if abs(term) > maxterm:
            maxterm = abs(term)
        if not math.isfinite(term):
            overflow = True
    result = 1.0 - total
    # only trust the float answer if it is finite and not swamped by cancellation
    if (not overflow) and math.isfinite(result) and result > maxterm * 1e-8:
        return result
    return coalesce_all_prob_mp(m, u)


## same probability but in arbitrary precision (fallback for the hard cases)
def coalesce_all_prob_mp(m, u):
    mp.dps = 50 + m   # precision grows with the number of lineages
    uu = mpf(u)
    rates = [mpf(j * (j - 1)) / 2 for j in range(2, m + 1)]
    total = mpf(0)
    for j in range(len(rates)):
        term = mpf(1)
        for i in range(len(rates)):
            if i != j:
                term = term * rates[i] / (rates[i] - rates[j])
        total = total + term * mpexp(-rates[j] * uu)
    return float(1 - total)


# bounded coalescent for branches feeding a duplication
## bounded coalescence so that duplication child lineages can merge before
## the duplication point, without us ending up with loads of polytomy
## from forcing them to all merge exactly at the duplication point
def bounded_next_time(k, T, rng):
    lam = k * (k - 1) / 2.0
    gT = coalesce_all_prob(k, T)
    if gT <= 0.0:   # safety net for extreme underflow (very rare)
        return T * rng.random()
    target = rng.random()
    lo, hi = 0.0, T
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        S = math.exp(-lam * mid) * coalesce_all_prob(k, T - mid) / gT
        if S > target:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


## run the bounded coalescent race, so that all merge and we dont get polytomy
def bounded_race(pool, start_time, top_time, theta, rng):
    time = start_time
    while len(pool) > 1:
        k = len(pool)
        wait = bounded_next_time(k, top_time - time, rng)
        time = time + wait
        merge_two(pool, time, theta, rng)
    return pool


# reading SaGePhy annotations
def node_id(name, tag):
    return name[len(tag):].split("_")[0]


################################################
# main entry point
################################################

def main(raw_tree, dup_ids, theta=1000, seed=None, tag="XTAGX"):
    """
    Rearrange raw_tree (an ete4 Tree, e.g. the parsed .pruned.tree) under the
    MSC at the given theta, treating dup_ids (a set of SaGePhy guest vertex
    IDs, e.g. {"6", "9"}) as fixed duplication points.

    Returns (coalescent_tree, norm_rf), where coalescent_tree is an ete4 Tree
    and norm_rf is the normalised Robinson-Foulds distance between raw_tree
    and coalescent_tree (0 = identical topology, i.e. no realised ILS).
    """
    if seed is None:
        seed = random.randint(0, 2**32 - 1)
    gene_tree = rearrange(raw_tree, theta, seed, dup_ids, tag)
    coalescent_tree = Tree(gene_tree, parser=1)
    rf = raw_tree.compare(coalescent_tree)["norm_rf"]
    return coalescent_tree, rf
