'''
Integration test for the first stage of --Get-Parameters --tool ORTHOXML:
turning an OrthoXML file plus proteomes into one FASTA file per orthogroup
and a species tree with numeric species names.
'''

import contextlib
import io
import os
import sys
import tempfile

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)

from Bio import SeqIO

import scripts.GetParameters.OrthoXML_MakeOrthogroupSequenceFiles as ox

XML = '''<?xml version="1.0"?>
<orthoXML xmlns="http://orthoXML.org/2011/" version="0.3" origin="t">
  <species name="Alpha" NCBITaxId="1"><database name="d" version="1"><genes>
    <gene id="1" protId="tr|A1|X"/><gene id="2" protId="A2"/><gene id="3" protId="A3"/>
  </genes></database></species>
  <species name="Beta" NCBITaxId="2"><database name="d" version="1"><genes>
    <gene id="4" protId="B1"/><gene id="5" protId="B2"/><gene id="6" protId="B9_missing"/>
  </genes></database></species>
  <groups>
    <orthologGroup id="HOG:1">
      <paralogGroup><geneRef id="1"/><geneRef id="2"/></paralogGroup>
      <geneRef id="4"/>
    </orthologGroup>
    <orthologGroup id="HOG:2"><geneRef id="3"/><geneRef id="5"/></orthologGroup>
    <orthologGroup id="HOG:3"><geneRef id="6"/></orthologGroup>
  </groups>
</orthoXML>
'''

ALPHA = ">tr|A1|X description\nMKVLAA\n>A2\nMKVLGG\n>A3\nWWWWWW\n"
BETA = ">B1\nMKVLTT\n>B2\nYYYYYY\n"


def build(tmp):
    prot = os.path.join(tmp, "prot")
    os.makedirs(prot)
    for name, text in (("Alpha.fa", ALPHA), ("Beta.fa", BETA)):
        with open(os.path.join(prot, name), "w") as f:
            f.write(text)
    xml = os.path.join(tmp, "x.orthoxml")
    with open(xml, "w") as f:
        f.write(XML)
    tree = os.path.join(tmp, "tree.nwk")
    with open(tree, "w") as f:
        f.write("(Alpha:1,Beta:1);")
    return xml, prot, tree


def run_main(tmp):
    xml, prot, tree = build(tmp)
    out = os.path.join(tmp, "out")
    os.makedirs(out)
    log = io.StringIO()
    with contextlib.redirect_stdout(log):
        ox.main(xml, prot, out, 1, tree)
    return out, log.getvalue()


def test_every_gene_of_a_group_is_written_with_its_species_number():
    with tempfile.TemporaryDirectory() as tmp:
        out, _ = run_main(tmp)
        seqs = {r.id: str(r.seq) for r in SeqIO.parse(
            os.path.join(out, "TrainingResults", "Orthogroup_Sequences", "OG_0000001.fa"), "fasta")}
    # Alpha is species 1 and Beta is species 2; the "|" in the accession is made safe
    assert seqs == {"1_tr_A1_X": "MKVLAA", "1_A2": "MKVLGG", "2_B1": "MKVLTT"}


def test_both_paralogs_are_kept():
    with tempfile.TemporaryDirectory() as tmp:
        out, _ = run_main(tmp)
        path = os.path.join(out, "TrainingResults", "Orthogroup_Sequences", "OG_0000001.fa")
        ids = [r.id for r in SeqIO.parse(path, "fasta")]
    assert sum(i.startswith("1_") for i in ids) == 2


def test_groups_with_fewer_than_two_resolved_genes_are_skipped():
    with tempfile.TemporaryDirectory() as tmp:
        out, log = run_main(tmp)
        files = sorted(os.listdir(os.path.join(out, "TrainingResults", "Orthogroup_Sequences")))
    assert files == ["OG_0000001.fa", "OG_0000002.fa"]
    assert "1 OrthoXML gene entries were not found" in log        # gene 6 has no sequence


def test_species_tree_is_copied_and_renumbered():
    with tempfile.TemporaryDirectory() as tmp:
        out, _ = run_main(tmp)
        numeric = open(os.path.join(out, "Species_Tree", "SpeciesTree_rooted_node_labels_numeric.txt")).read()
        mapping = open(os.path.join(out, "Species_Tree", "species_name_to_number.csv")).read().split()
    assert "1:1" in numeric.replace(" ", "") and "2:1" in numeric.replace(" ", "")
    assert mapping == ["species,number", "Alpha,1", "Beta,2"]


def test_a_species_in_the_orthoxml_without_a_proteome_is_an_error():
    with tempfile.TemporaryDirectory() as tmp:
        xml, prot, tree = build(tmp)
        os.remove(os.path.join(prot, "Beta.fa"))
        out = os.path.join(tmp, "out")
        os.makedirs(out)
        try:
            with contextlib.redirect_stdout(io.StringIO()):
                ox.main(xml, prot, out, 1, tree)
        except ValueError as err:
            assert "Beta" in str(err)
            return
    raise AssertionError("a missing proteome was accepted")


def test_a_species_tree_is_required():
    with tempfile.TemporaryDirectory() as tmp:
        xml, prot, tree = build(tmp)
        try:
            with contextlib.redirect_stdout(io.StringIO()):
                ox.main(xml, prot, os.path.join(tmp, "out"), 1, None)
        except ValueError as err:
            assert "tool-tree" in str(err)
            return
    raise AssertionError("a missing species tree was accepted")
