#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
OrthoXML GetParameters sequence-file builder.

Generic input adapter: consumes any orthology result exported in the
OrthoXML standard (https://orthoxml.org) rather than a tool-specific output
format, so tools that natively emit OrthoXML (or that can be converted to
it, e.g. ProteinOrtho via proteinortho2xml.pl) can plug into GetParameters
without a bespoke parser each.

Deliberately narrow in scope: this script only uses OrthoXML for
root-level orthogroup ("root HOG") membership. It does NOT use the
hierarchical HOG structure for anything else (no `to-nhx`, no HOG-native
duplication/rooting) - gene-ID conventions and tree quality are not
guaranteed consistent across OrthoXML producers (the spec itself makes
geneId/protId/transcriptId on <gene> all optional, and different tools
populate different subsets in different formats). Everything downstream
of this script - alignment, tree estimation, duplication calling via
CallDuplications.py's species-tree-guided rerooting - runs exactly as it
does for every other tool in this pipeline.

IMPORTANT: root-HOG membership is extracted with our own streaming parser
(parse_orthoxml, below), NOT via orthoxml-tools' `export-ogs` CLI.
export-ogs computes a "maximum orthologous group": at every <paralogGroup>
duplication node it keeps only the single branch with the most distinct
species and silently discards the rest, and for a leaf-level
<paralogGroup> (direct duplicate geneRefs with no further nesting) it
keeps only the FIRST geneRef. That is the wrong data source for this
pipeline - we need every paralog/in-paralog preserved so alignment,
tree-building and CallDuplications.py's own rerooting can rediscover
duplications themselves, the same way they do for every other tool's
input. Confirmed against hand-built OrthoXML files with (a) a
leaf-level in-species duplication and (b) a duplication predating
speciation: export-ogs silently drops 1-2 genes per case, which is the
root cause of duplication counts coming back as zero downstream - the
sequence files it would have driven this script from had already had
their duplicated genes pruned out before a gene tree was ever built.

Species assignment is likewise read directly from the OrthoXML file's own
<species name="..."><genes><gene id="..."/></genes></species> structure
(authoritative - this is what the tool that produced the file actually
asserts), rather than being inferred sideways from which proteome FASTA
file happens to contain a sequence whose header matches a gene's id.
The proteome FASTA files are used ONLY to fetch the actual sequence for
a gene once its species is already known from the OrthoXML file; sequence
lookup is scoped to that one species' own file, trying each of protId,
geneId, id (in that order) as an exact match against FASTA headers.
"""

import os
import csv
import shutil
from lxml import etree

from Bio import SeqIO
from ete4 import Tree


# candidate gene-ID fields to try, in order, when matching an OrthoXML gene
# entry against a proteome FASTA header. "id" is always present (it's the
# internal reference OrthoXML uses for geneRef linking); protId/geneId are
# optional per the spec and producer-dependent.
ID_FIELDS = ("protId", "geneId", "id")


def _local_tag(tag):
    """Strip the XML namespace off an lxml element tag, e.g.
    '{http://orthoXML.org/2011/}species' -> 'species'."""
    return tag.split("}", 1)[-1]


def _collect_generefs(elem):
    """
    Recursively collect every <geneRef id="..."/> nested anywhere under
    elem - a full, unpruned union across all child <orthologGroup> and
    <paralogGroup> elements. Deliberately does NOT do any of export-ogs'
    "maximum orthologous group" branch-picking/truncation (see module
    docstring) - every paralog and in-paralog is kept.
    """
    ids = []
    for child in elem:
        tag = _local_tag(child.tag)
        if tag == "geneRef":
            gid = child.get("id")
            if gid is not None:
                ids.append(gid)
        elif tag in ("orthologGroup", "paralogGroup"):
            ids.extend(_collect_generefs(child))
    return ids


def parse_orthoxml(orthoxml_file):
    """
    Single streaming pass over the OrthoXML file (memory-safe for large
    files via incremental element clearing, same approach orthoxml-tools
    itself uses internally).

    Returns:
        orthogroups: {og_name: [gene_id, ...]}
            Full, unpruned membership of every top-level <orthologGroup>/
            <paralogGroup> ("root HOG"). og_name follows the same
            "OG_0000001"-style convention export-ogs used, for continuity.

        gene_attrs: {gene_id: {"species": species_name, "id":.., "protId":.., "geneId":..}}
            Built directly from each <species> block's own <gene> children
            - this is OrthoXML's own authoritative species assignment.
    """
    orthogroups = {}
    gene_attrs = {}

    og_depth = 0
    pg_depth = 0
    og_counter = 0

    context = etree.iterparse(orthoxml_file, events=("start", "end"))
    for event, elem in context:
        tag = _local_tag(elem.tag)

        if event == "start":
            if tag == "orthologGroup":
                og_depth += 1
            elif tag == "paralogGroup":
                pg_depth += 1
            continue

        # event == "end"
        if tag == "species":
            species_name = elem.get("name")
            for descendant in elem.iter():
                if _local_tag(descendant.tag) != "gene":
                    continue
                gid = descendant.get("id")
                if gid is None:
                    continue
                gene_attrs[gid] = {
                    "species": species_name,
                    "id": descendant.get("id"),
                    "protId": descendant.get("protId"),
                    "geneId": descendant.get("geneId"),
                }
            elem.clear()

        elif tag == "orthologGroup":
            og_depth -= 1
            if og_depth == 0 and pg_depth == 0:
                og_counter += 1
                orthogroups[f"OG_{og_counter:07d}"] = _collect_generefs(elem)
                elem.clear()

        elif tag == "paralogGroup":
            pg_depth -= 1
            if og_depth == 0 and pg_depth == 0:
                og_counter += 1
                orthogroups[f"OG_{og_counter:07d}"] = _collect_generefs(elem)
                elem.clear()

    return orthogroups, gene_attrs


def sanitise_gene_id(gene_id):
    """Make a gene accession safe to embed in a numeric FASTA/tree leaf name.

    Several downstream tools in this pipeline (aligners, tree builders,
    Newick parsers, CallDuplications.py's own leaf-name parsing) do not
    reliably preserve arbitrary characters in taxon/sequence names. Real
    UniProt/QfO-style proteome headers commonly look like
    "tr|B8BSL4|B8BSL4_THAPS" - the "|" characters are a known landmine for
    that toolchain (silently stripped/mangled by one or more steps between
    alignment and tree building), which breaks the exact-string match
    CallDuplications.py relies on to map gene-tree leaves back to species,
    and duplication counts come back silently empty rather than raising an
    error. Kept as a separate transform from the raw accession (which is
    still used, unmodified, as the gene_index lookup key against the
    original FASTA headers) so exact-match resolution is unaffected -
    sanitisation only applies to the ID written into the numeric output
    FASTA/tree files."""
    safe = "".join(
        ch if (ch.isalnum() or ch in ("_", "-")) else "_"
        for ch in gene_id
    )
    return safe or "GENE_UNKNOWN"


def sanitise_og_name(og_name):
    """Make an orthogroup name safe to use as a filename, without assuming
    any particular source naming convention (e.g. FastOMA's "OG_0000001" vs
    some other producer's own scheme)."""
    safe = "".join(
        ch if (ch.isalnum() or ch in ("_", "-")) else "_"
        for ch in og_name
    )
    return safe or "OG_UNKNOWN"


def index_proteomes(proteome_dir):
    """
    Index every proteome FASTA in proteome_dir, assigning a numeric
    species ID per file (sorted filename order - same convention used by
    the other *_MakeOrthogroupSequenceFiles.py scripts in this repo).

    Unlike the earlier version of this script, the sequence index is kept
    PER SPECIES rather than as one combined cross-species dict. Species is
    now known independently from the OrthoXML file itself (parse_orthoxml's
    gene_attrs), so there is no need to assume accessions are globally
    unique across every species in the dataset - only unique within a
    single species' own proteome file, which is a much weaker and safer
    assumption. A duplicate accession appearing in two different species'
    files is completely fine now; it simply resolves independently within
    whichever species OrthoXML says the gene belongs to.

    Returns:
        species_to_id: {species_name: numeric_id}
        gene_index: {species_name: {accession: SeqRecord}}
    """
    species_to_id = {}
    gene_index = {}

    proteome_files = sorted(os.listdir(proteome_dir))
    species_counter = 0
    total_sequences = 0

    for proteome_file in proteome_files:
        if not proteome_file.endswith((".fa", ".faa", ".fasta")):
            continue

        species_name = os.path.splitext(proteome_file)[0]
        species_counter += 1
        species_to_id[species_name] = species_counter

        file_path = os.path.join(proteome_dir, proteome_file)
        per_species_index = gene_index.setdefault(species_name, {})

        for record in SeqIO.parse(file_path, "fasta"):
            accession = record.id  # Bio.SeqIO already splits off the description

            if accession in per_species_index:
                raise ValueError(
                    f"Duplicate gene accession '{accession}' found twice "
                    f"within proteome '{species_name}' ({file_path}). "
                    "Sequence resolution assumes accessions are unique "
                    "within a single species' proteome file."
                )

            per_species_index[accession] = record
            total_sequences += 1

    if total_sequences == 0:
        entries = sorted(os.listdir(proteome_dir))
        subdirs = [
            e for e in entries if os.path.isdir(os.path.join(proteome_dir, e))
        ]
        hint = ""
        if subdirs:
            hint = (
                f" --tools-proteomes must point directly at a flat directory "
                f"of .fa/.faa/.fasta files - it does not search subdirectories. "
                f"Found subdirector{'y' if len(subdirs) == 1 else 'ies'} instead: "
                f"{subdirs[:10]}. Some tools (e.g. FastOMA) nest proteomes one "
                f"level down (such as <input_folder>/proteome/) - check whether "
                f"one of these subdirectories is what you meant to pass."
            )
        raise ValueError(
            f"No protein sequences were indexed from: {proteome_dir} "
            f"(found {len(entries)} entries, none ending in .fa/.faa/.fasta "
            f"directly in this directory).{hint}"
        )

    if not species_to_id:
        raise ValueError(f"No species FASTA files were found in: {proteome_dir}")

    return species_to_id, gene_index


def resolve_gene(gene_id, gene_attrs, gene_index):
    """
    Look up the species OrthoXML asserts for gene_id, then try each
    candidate identifier (protId, geneId, id - in ID_FIELDS order) as an
    exact match against that ONE species' own proteome index. Species
    assignment is always taken from OrthoXML regardless of which field
    (if any) ends up resolving the sequence - the two are decoupled.

    Returns:
        (species_name, SeqRecord) on success, or (None, reason_str) where
        reason_str explains why resolution failed (for warning messages).
    """
    attrs = gene_attrs.get(gene_id)
    if attrs is None:
        return None, "gene id from <groups> was never declared in any <species> block"

    species_name = attrs["species"]
    per_species_index = gene_index.get(species_name)
    if per_species_index is None:
        return None, f"species '{species_name}' has no matching proteome file"

    for field in ID_FIELDS:
        value = attrs.get(field)
        if not value:
            continue
        if value in per_species_index:
            return species_name, per_species_index[value]

    return None, (
        f"none of protId/geneId/id matched a sequence in species "
        f"'{species_name}''s proteome file"
    )


def save_numeric_species_tree(
    species_tree_src,
    numeric_species_tree_path,
    species_mapping_path,
    species_to_id,
):
    """
    Save a copy of the species tree with leaf names replaced by the same
    numeric species IDs used in the orthogroup sequence FASTA files.
    """
    tree = Tree(species_tree_src, parser=1)

    for leaf in tree.leaves():
        if leaf.name not in species_to_id:
            raise ValueError(
                f"Species '{leaf.name}' is in the species tree but was not "
                "found in the proteome directory."
            )
        leaf.name = str(species_to_id[leaf.name])

    tree.write(outfile=numeric_species_tree_path, parser=1)

    with open(species_mapping_path, "w", newline="") as out:
        writer = csv.writer(out)
        writer.writerow(["species", "number"])
        for species in sorted(species_to_id):
            writer.writerow([species, species_to_id[species]])

    print(f"Wrote numeric species tree to {numeric_species_tree_path}")


def main(orthoxml_file, proteome_dir, output_base, n_threads, tool_tree):
    """
    Parameters
    ----------
    orthoxml_file : str
        Path to the .orthoxml file to build orthogroup sequence files from.

    proteome_dir : str
        Directory of per-species proteome FASTA files (passed by
        OrthoSimulator.py via --tools-proteomes).

    output_base : str
        OrthoSimulator output folder.

    n_threads : int
        Accepted for interface consistency with the other Train scripts;
        this step's own work (writing per-OG FASTA files) is not itself
        parallelised, matching FastOMA_MakeOrthogroupSequenceFiles.py.

    tool_tree : str
        Path to the species tree (required - OrthoXML files are not
        guaranteed to carry a species tree with usable branch lengths,
        same requirement as Broccoli/SonicParanoid2).
    """
    orthoxml_file = os.path.abspath(orthoxml_file)
    proteome_dir = os.path.abspath(proteome_dir)
    output_base = os.path.abspath(output_base)

    if not os.path.isfile(orthoxml_file):
        raise ValueError(f"OrthoXML file not found: {orthoxml_file}")
    if not os.path.isdir(proteome_dir):
        raise ValueError(f"Proteome directory not found: {proteome_dir}")
    if not tool_tree:
        raise ValueError(
            "--tool-tree is required for ORTHOXML (OrthoXML files are not "
            "guaranteed to carry a species tree with usable branch lengths)."
        )
    tool_tree = os.path.abspath(tool_tree)
    if not os.path.isfile(tool_tree):
        raise ValueError(f"Species tree not found: {tool_tree}")

    seq_dir = os.path.join(output_base, "TrainingResults", "Orthogroup_Sequences")
    # NOTE: unlike the other *_MakeOrthogroupSequenceFiles.py scripts, which
    # write Species_Tree/ into the raw tool's own output directory, ORTHOXML's
    # "tool output" is a single file with no directory of its own to write
    # into (and that directory may not even be writable). Species_Tree/ is
    # written under output_base instead, which OrthoSimulator always creates and
    # owns. OrthoXMLTrain.py passes output_base as CallDuplications.py's
    # "tool output directory" argument to match.
    species_tree_dir = os.path.join(output_base, "Species_Tree")

    os.makedirs(seq_dir, exist_ok=True)
    os.makedirs(species_tree_dir, exist_ok=True)

    species_tree_dst = os.path.join(species_tree_dir, "SpeciesTree_rooted_node_labels.txt")
    numeric_species_tree_path = os.path.join(
        species_tree_dir, "SpeciesTree_rooted_node_labels_numeric.txt"
    )
    species_mapping_path = os.path.join(species_tree_dir, "species_name_to_number.csv")

    print(f"OrthoXML file: {orthoxml_file}")
    print(f"Proteomes: {proteome_dir}")
    print(f"Sequence output: {seq_dir}")
    print(f"Species tree: {tool_tree}")

    # ============================================================
    # Root-orthogroup membership + authoritative species assignment,
    # both read directly from the OrthoXML file (see module docstring for
    # why export-ogs is NOT used here - it prunes duplications)
    # ============================================================

    orthogroups, gene_attrs = parse_orthoxml(orthoxml_file)
    total_genes_in_groups = sum(len(ids) for ids in orthogroups.values())
    print(
        f"Loaded {len(orthogroups)} root orthogroups from OrthoXML "
        f"({total_genes_in_groups} gene references, unpruned)"
    )

    # ============================================================
    # Index proteome sequences and assign numeric species IDs
    # ============================================================

    species_to_id, gene_index = index_proteomes(proteome_dir)
    total_sequences = sum(len(idx) for idx in gene_index.values())
    print(f"Indexed {total_sequences} sequences across {len(species_to_id)} species")

    # Cross-check: every species OrthoXML declares must have a matching
    # proteome file, or sequence resolution can never succeed for its genes.
    orthoxml_species = {attrs["species"] for attrs in gene_attrs.values()}
    proteome_species = set(species_to_id)
    missing_species = sorted(orthoxml_species - proteome_species)
    if missing_species:
        raise ValueError(
            f"{len(missing_species)} species named in the OrthoXML file have "
            f"no matching proteome FASTA file in {proteome_dir} (species name "
            "must exactly match the FASTA filename minus extension). "
            f"Missing (up to 10 shown): {missing_species[:10]}. "
            f"Proteome species found (up to 10 shown): {sorted(proteome_species)[:10]}."
        )

    # ============================================================
    # Write orthogroup sequence files
    # ============================================================

    written = 0
    missing_genes = 0
    unresolved_examples = []
    # Tracks sanitised numeric IDs across the WHOLE run (not just within one
    # orthogroup) - sanitisation collapses any non-alphanumeric character to
    # "_", which could in principle map two distinct accessions to the same
    # string. Vanishingly unlikely for real accessions, but silently letting
    # two different genes collide would corrupt output, so it's checked.
    used_numeric_ids = {}

    for og_name, gene_ids in orthogroups.items():
        og_safe = sanitise_og_name(og_name)

        # dedupe on the OrthoXML gene id itself (always unique per <gene>
        # element) - a gene id can in principle appear more than once
        # inside one root HOG's tree structure. Deliberately NOT deduping
        # on the resolved accession string: accessions are only guaranteed
        # unique within a single species now (see index_proteomes), so the
        # same accession string can legitimately appear in two different
        # species within one orthogroup and both must be kept.
        seen_gene_ids = set()
        records_num = []

        for gene_id in gene_ids:
            if gene_id in seen_gene_ids:
                continue
            seen_gene_ids.add(gene_id)

            species_name, resolved = resolve_gene(gene_id, gene_attrs, gene_index)
            if species_name is None:
                missing_genes += 1
                if len(unresolved_examples) < 10:
                    unresolved_examples.append((og_name, gene_id, resolved))
                continue

            record = resolved
            species_id = species_to_id[species_name]
            numeric_id = f"{species_id}_{sanitise_gene_id(record.id)}"
            if numeric_id in used_numeric_ids and used_numeric_ids[numeric_id] != record.id:
                raise ValueError(
                    f"Sanitised gene ID collision: '{numeric_id}' was produced "
                    f"for two different original accessions "
                    f"('{used_numeric_ids[numeric_id]}' and '{record.id}'). "
                    "Gene-ID sanitisation replaces every non-alphanumeric "
                    "character with '_', which can in principle collide for "
                    "unusual accessions - this dataset apparently triggers "
                    "that. A different sanitisation scheme would be needed."
                )
            used_numeric_ids[numeric_id] = record.id

            rec_num = record[:]
            rec_num.id = numeric_id
            rec_num.name = rec_num.id
            rec_num.description = ""
            records_num.append(rec_num)

        if len(records_num) < 2:
            continue

        out_fasta = os.path.join(seq_dir, f"{og_safe}.fa")
        SeqIO.write(records_num, out_fasta, "fasta")
        written += 1

    print(f"Wrote {written} orthogroup sequence files")

    if missing_genes > 0:
        print(
            f"WARNING: {missing_genes} OrthoXML gene entries were not found "
            "in the indexed proteomes (tried protId, geneId, id, in that "
            "order, exact match only)."
        )
        for og_name, gene_id, reason in unresolved_examples:
            print(f"  example: group={og_name} gene_id={gene_id} reason={reason}")

    if written == 0:
        raise ValueError(
            "No orthogroup sequence files were written. Check that the "
            "OrthoXML gene identifiers (protId/geneId/id) match the "
            "proteome FASTA headers."
        )

    # ============================================================
    # Species tree
    # ============================================================

    shutil.copyfile(tool_tree, species_tree_dst)
    print(f"Copied species tree to {species_tree_dst}")

    save_numeric_species_tree(
        species_tree_src=tool_tree,
        numeric_species_tree_path=numeric_species_tree_path,
        species_mapping_path=species_mapping_path,
        species_to_id=species_to_id,
    )

    print("Sequence generation complete.")
