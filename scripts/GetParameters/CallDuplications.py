#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
from copy import deepcopy
from multiprocessing import Pool
#from ete3 import Tree
from ete4 import Tree


# ============================================================
# Helpers
# ============================================================

def GetSpeciesID(gene_name):
    """
    Gene names are assumed to look like:
        1_geneA
        12_geneB

    Species ID is the part before the first underscore.
    """
    return gene_name.split("_")[0]


def GetRepresentedSpecies(gene_tree):
    species = set()

    for leaf in gene_tree.leaves():
        species.add(GetSpeciesID(leaf.name))

    return species


def GetGenesFromSpecies(gene_tree, species_set):
    genes = []

    for leaf in gene_tree.leaves():
        species_id = GetSpeciesID(leaf.name)

        if species_id in species_set:
            genes.append(leaf.name)

    return genes


def GetSpeciesUnderGeneNode(gene_node):
    species = set()

    for leaf in gene_node.leaves():
        species.add(GetSpeciesID(leaf.name))

    return species


def GetSpeciesUnderSpeciesNode(species_node):
    #print(species_node.leaves())
    #print("species node")
    #print(species_node.write(parser=1))
    #print(species_node.leaf_names())
    return set(species_node.leaf_names())


def JoinSet(values):
    """
    Write sets/lists cleanly in troubleshooting output.
    """
    if not values:
        return ""

    return ",".join(sorted(str(v) for v in values))


# ============================================================
# Rooting: OrthoFinder-style species-tree-guided rooting
# ============================================================

def FindRelevantSpeciesSplits(species_tree, represented_species):
    """
    Find the first relevant split in the rooted species tree.

    Start at the species-tree root.
    Move down while all represented species are in only one child.
    Stop when represented species occur in at least two child clades.

    For a binary node:
        A | B

    For a polytomy:
        A | rest
        B | rest
        C | rest
    """
    node = species_tree
    

    while True:
        children = node.get_children()
        #print("children")
        #print(children)
        
        if not children:
            return []

        child_species_sets = [
            GetSpeciesUnderSpeciesNode(child)
            for child in children
        ]
        #print(child_species_sets)
        #print("###")
        #print(represented_species)
        have = [
            len(child_species & represented_species) > 0
            for child_species in child_species_sets
        ]
        #print(have)
        if sum(have) >= 2:
            break

        if sum(have) == 0:
            return []

        node = children[have.index(True)]

    splits = []

    #print(child_species_sets)
    #print("###")
    for child_species in child_species_sets:
        outgroup_species = child_species & represented_species
        ingroup_species = represented_species - outgroup_species

        if not outgroup_species:
            continue

        if not ingroup_species:
            continue
        
        splits.append({
            "species_lca_name": node.name,
            "outgroup_species": outgroup_species,
            "ingroup_species": ingroup_species,
        })

    return splits


def LabelSpeciesSet(species_set, outgroup_species, ingroup_species):
    """
    Label a species set relative to a candidate outgroup/ingroup split.

    OUT:
        contains outgroup species only

    IN:
        contains ingroup species only

    MIXED:
        contains both outgroup and ingroup species

    EMPTY:
        contains neither
    """
    has_out = len(species_set & outgroup_species) > 0
    has_in = len(species_set & ingroup_species) > 0

    if has_out and has_in:
        return "MIXED"

    if has_out:
        return "OUT"

    if has_in:
        return "IN"

    return "EMPTY"


def ScoreRootCandidate(
    down_species,
    up_species,
    outgroup_species,
    ingroup_species,
):
    """
    Score an imperfect root candidate.

    Higher is better.

    Perfect clean split:
        one side is OUT
        the other side is IN

    If no clean split exists, this rewards branches that best separate
    outgroup and ingroup species.
    """
    down_out = len(down_species & outgroup_species)
    down_in = len(down_species & ingroup_species)

    up_out = len(up_species & outgroup_species)
    up_in = len(up_species & ingroup_species)

    n_out = max(len(outgroup_species), 1)
    n_in = max(len(ingroup_species), 1)

    # Orientation 1:
    # down = outgroup, up = ingroup
    score1 = (down_out / n_out) * (up_in / n_in)

    # Orientation 2:
    # down = ingroup, up = outgroup
    score2 = (down_in / n_in) * (up_out / n_out)

    return max(score1, score2)


def FindBestGeneRootForSplit(gene_tree, represented_species, split):
    """
    Find the best root point in the gene tree for one species-tree split.

    For every non-root node, compare:
        species below node
        species above node

    If one side is OUT and the other is IN, this is a clean root.

    If no clean root exists, keep the highest-scoring imperfect root.
    """
    outgroup_species = split["outgroup_species"]
    ingroup_species = split["ingroup_species"]

    best_node = None
    best_score = -1
    best_down_label = ""
    best_up_label = ""
    best_down_species = set()
    best_up_species = set()

    for node in gene_tree.traverse("postorder"):
        if node.is_root:
            continue

        down_species = GetSpeciesUnderGeneNode(node)
        up_species = represented_species - down_species

        down_label = LabelSpeciesSet(
            down_species,
            outgroup_species,
            ingroup_species,
        )

        up_label = LabelSpeciesSet(
            up_species,
            outgroup_species,
            ingroup_species,
        )

        clean = (
            (down_label == "OUT" and up_label == "IN") or
            (down_label == "IN" and up_label == "OUT")
        )

        if clean:
            return {
                "node": node,
                "is_clean": True,
                "score": 1.0,
                "down_label": down_label,
                "up_label": up_label,
                "down_species": down_species,
                "up_species": up_species,
                "outgroup_species": outgroup_species,
                "ingroup_species": ingroup_species,
                "species_lca_name": split["species_lca_name"],
            }

        score = ScoreRootCandidate(
            down_species,
            up_species,
            outgroup_species,
            ingroup_species,
        )

        if score > best_score:
            best_node = node
            best_score = score
            best_down_label = down_label
            best_up_label = up_label
            best_down_species = down_species
            best_up_species = up_species

    if best_node is None:
        return None

    return {
        "node": best_node,
        "is_clean": False,
        "score": best_score,
        "down_label": best_down_label,
        "up_label": best_up_label,
        "down_species": best_down_species,
        "up_species": best_up_species,
        "outgroup_species": outgroup_species,
        "ingroup_species": ingroup_species,
        "species_lca_name": split["species_lca_name"],
    }


def FindBestGeneRoot(gene_tree, species_tree, represented_species):
    """
    Try all relevant species-tree child-vs-rest splits.

    Prefer the first clean gene-tree split.

    If no clean split exists, choose the highest-scoring imperfect split.
    """
    splits = FindRelevantSpeciesSplits(
        species_tree,
        represented_species,
    )
    #print(splits)
    if not splits:
        return None, 0

    best = None

    for split in splits:
        candidate = FindBestGeneRootForSplit(
            gene_tree,
            represented_species,
            split,
        )

        if candidate is None:
            continue

        if candidate["is_clean"]:
            return candidate, len(splits)

        if best is None or candidate["score"] > best["score"]:
            best = candidate

    return best, len(splits)


def RootGeneTreeWithSpeciesTree(gene_tree, species_tree):
    """
    Root one gene tree using an OrthoFinder-style species-tree-guided method.

    Returns:
        rooted_tree, info
    """
    rooted_tree = deepcopy(gene_tree)

    represented_species = GetRepresentedSpecies(rooted_tree)
    species_tree_species = set(species_tree.leaf_names())

    info = {
        "status": "",
        "represented_species": represented_species,
        "missing_species": set(),
        "species_lca_name": "",
        "outgroup_species": set(),
        "ingroup_species": set(),
        "outgroup_genes": [],
        "outgroup_node_name": "",
        "rooting_method": "",
        "candidate_count": 0,
        "is_monophyletic": "",
        "extra_species": set(),
        "down_label": "",
        "up_label": "",
        "down_species": set(),
        "up_species": set(),
        "root_score": "",
        "error": "",
    }

    missing_species = represented_species - species_tree_species
    info["missing_species"] = missing_species

    if missing_species:
        info["status"] = "unrootable_species_missing_from_species_tree"
        return rooted_tree, info

    if len(represented_species) < 2:
        info["status"] = "unrootable_single_species"
        return rooted_tree, info

    root_candidate, n_splits = FindBestGeneRoot(
        rooted_tree,
        species_tree,
        represented_species,
    )
    info["candidate_count"] = n_splits

    if root_candidate is None:
        info["status"] = "unrootable_no_gene_root_candidate"
        return rooted_tree, info

    outgroup_node = root_candidate["node"]

    try:
        rooted_tree.set_outgroup(outgroup_node)

    except Exception as err:
        info["status"] = "rooting_failed"
        info["error"] = str(err)
        return rooted_tree, info

    info["species_lca_name"] = root_candidate["species_lca_name"]
    info["outgroup_species"] = root_candidate["outgroup_species"]
    info["ingroup_species"] = root_candidate["ingroup_species"]

    info["outgroup_genes"] = GetGenesFromSpecies(
        rooted_tree,
        root_candidate["outgroup_species"],
    )

    info["outgroup_node_name"] = outgroup_node.name
    info["down_label"] = root_candidate["down_label"]
    info["up_label"] = root_candidate["up_label"]
    info["down_species"] = root_candidate["down_species"]
    info["up_species"] = root_candidate["up_species"]
    info["root_score"] = root_candidate["score"]

    if root_candidate["is_clean"]:
        info["status"] = "rooted_clean_split"
        info["rooting_method"] = "orthofinder_style_clean_split"
        info["is_monophyletic"] = "True"
    else:
        info["status"] = "rooted_best_scoring_split"
        info["rooting_method"] = "orthofinder_style_best_score"
        info["is_monophyletic"] = "False"

    return rooted_tree, info


# ============================================================
# Duplication counting
# ============================================================

def CountSpeciesOverlapDuplications(rooted_gene_tree):
    """
    Count duplications using the species-overlap method.

    One internal node = one duplication if any species appears in
    two or more child clades.
    """
    duplication_count = 0

    for node in rooted_gene_tree.traverse("postorder"):
        if node.is_leaf:
            continue

        child_species_sets = []

        for child in node.children:
            species_set = set()

            for leaf in child.leaves():
                species_set.add(GetSpeciesID(leaf.name))

            child_species_sets.append(species_set)

        is_duplication = False

        for i in range(len(child_species_sets)):
            for j in range(i + 1, len(child_species_sets)):
                if child_species_sets[i] & child_species_sets[j]:
                    is_duplication = True
                    break

            if is_duplication:
                break

        if is_duplication:
            duplication_count += 1

    return duplication_count


# ============================================================
# Worker
# ============================================================

def ProcessGeneTree(tree_file, species_tree_path):
    """
    Worker for multiprocessing.

    Returns a dictionary with both:
        - minimal duplication count info
        - troubleshooting info
    """
    orthogroup = os.path.basename(os.path.dirname(tree_file))

    try:
        species_tree = Tree(species_tree_path, parser=1)
        gene_tree = Tree(tree_file, parser=1)

        rooted_tree, info = RootGeneTreeWithSpeciesTree(
            gene_tree,
            species_tree,
        )

        duplication_count = CountSpeciesOverlapDuplications(rooted_tree)

        return {
            "orthogroup": orthogroup,
            "duplications": duplication_count,
            "status": info["status"],
            "represented_species": JoinSet(info["represented_species"]),
            "missing_species": JoinSet(info["missing_species"]),
            "species_lca_name": info["species_lca_name"],
            "outgroup_species": JoinSet(info["outgroup_species"]),
            "ingroup_species": JoinSet(info["ingroup_species"]),
            "outgroup_genes": JoinSet(info["outgroup_genes"]),
            "outgroup_node_name": info["outgroup_node_name"],
            "rooting_method": info["rooting_method"],
            "candidate_count": info["candidate_count"],
            "is_monophyletic": info["is_monophyletic"],
            "extra_species": JoinSet(info["extra_species"]),
            "down_label": info["down_label"],
            "up_label": info["up_label"],
            "down_species": JoinSet(info["down_species"]),
            "up_species": JoinSet(info["up_species"]),
            "root_score": info["root_score"],
            "error": info["error"],
        }

    except Exception as err:
        return {
            "orthogroup": orthogroup,
            "duplications": "NA",
            "status": "failed",
            "represented_species": "",
            "missing_species": "",
            "species_lca_name": "",
            "outgroup_species": "",
            "ingroup_species": "",
            "outgroup_genes": "",
            "outgroup_node_name": "",
            "rooting_method": "",
            "candidate_count": "",
            "is_monophyletic": "",
            "extra_species": "",
            "down_label": "",
            "up_label": "",
            "down_species": "",
            "up_species": "",
            "root_score": "",
            "error": str(err),
        }


# ============================================================
# Output
# ============================================================

def WriteDuplicationsFile(out_file, results):
    """
    Minimal useful file.

    Literally just one duplication count per line.
    """
    with open(out_file, "w") as out:
        for row in results:
            out.write(f"{row['duplications']}\n")


def WriteTroubleshootingFile(out_file, results):
    """
    Detailed troubleshooting file.
    """
    fields = [
        "orthogroup",
        "duplications",
        "status",
        "represented_species",
        "missing_species",
        "species_lca_name",
        "outgroup_species",
        "ingroup_species",
        "outgroup_genes",
        "outgroup_node_name",
        "rooting_method",
        "candidate_count",
        "is_monophyletic",
        "extra_species",
        "down_label",
        "up_label",
        "down_species",
        "up_species",
        "root_score",
        "error",
    ]

    with open(out_file, "w") as out:
        out.write("\t".join(fields) + "\n")

        for row in results:
            values = [
                str(row[field])
                for field in fields
            ]
            out.write("\t".join(values) + "\n")


# ============================================================
# Main callable function
# ============================================================

def main(fastoma_out, output_folder_path, n_threads):
    """
    Callable from the main pipeline.

    Inputs:
        fastoma_out:
            FastOMA output folder.

        output_folder_path:
            Main OrthoSimulator/output folder containing gene trees.

        n_threads:
            Number of processes.

    Expected species tree:
        fastoma_out/Species_Tree/SpeciesTree_rooted_node_labels_numeric.txt

    Expected gene-tree files:
        output_folder_path/TrainingResults/Orthogroup_analysis/*/*.treefile

    Outputs:
        output_folder_path/TrainingResults/Orthogroup_analysis/duplications.txt
        output_folder_path/TrainingResults/Orthogroup_analysis/duplication_rooting_debug.txt
    """

    species_tree_path = os.path.join(
        os.path.abspath(fastoma_out),
        "Species_Tree",
        "SpeciesTree_rooted_node_labels_numeric.txt",
    )

    gene_tree_dir = os.path.join(
        os.path.abspath(output_folder_path),
        "TrainingResults",
        "Orthogroup_analysis",
    )

    out_dir = gene_tree_dir

    duplications_file = os.path.join(
        out_dir,
        "duplications.txt",
    )

    debug_file = os.path.join(
        out_dir,
        "duplication_rooting_debug.txt",
    )

    if not os.path.isfile(species_tree_path):
        raise ValueError(f"Numeric species tree not found: {species_tree_path}")

    if not os.path.isdir(gene_tree_dir):
        raise ValueError(f"Gene tree directory not found: {gene_tree_dir}")

    os.makedirs(out_dir, exist_ok=True)

    tree_files = []

    for root, dirs, files in os.walk(gene_tree_dir):
        for name in files:
            if name.endswith(".treefile"):
                tree_files.append(os.path.join(root, name))

    tree_files = sorted(tree_files)

    print(f"Numeric species tree: {species_tree_path}")
    print(f"Gene tree folder:     {gene_tree_dir}")
    print(f"Duplication output:   {duplications_file}")
    print(f"Debug output:         {debug_file}")
    print(f"Gene trees found:     {len(tree_files)}")
    print(f"Threads:              {n_threads}")

    args = [
        (tree_file, species_tree_path)
        for tree_file in tree_files
    ]

    if int(n_threads) == 1:
        results = [
            ProcessGeneTree(*arg)
            for arg in args
        ]
    else:
        with Pool(int(n_threads)) as pool:
            results = pool.starmap(ProcessGeneTree, args)

    results = sorted(
        results,
        key=lambda row: row["orthogroup"],
    )

    WriteDuplicationsFile(duplications_file, results)
    WriteTroubleshootingFile(debug_file, results)

    print(f"Wrote duplication counts to {duplications_file}")
    print(f"Wrote troubleshooting info to {debug_file}")
