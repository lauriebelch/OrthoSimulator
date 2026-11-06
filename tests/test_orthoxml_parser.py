'''
Unit tests for the OrthoXML reader used by --Get-Parameters --tool ORTHOXML
(scripts/GetParameters/OrthoXML_MakeOrthogroupSequenceFiles.py).

The reader must keep every gene of a top-level group, including paralogs
and in-paralogs, and take each gene's species from the OrthoXML species block.
'''

import os
import sys
import tempfile

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)

import scripts.GetParameters.OrthoXML_MakeOrthogroupSequenceFiles as ox

NS = "http://orthoXML.org/2011/"

XML = f'''<?xml version="1.0" encoding="UTF-8"?>
<orthoXML xmlns="{NS}" version="0.3" origin="test">
  <species name="Homo_sapiens" NCBITaxId="9606">
    <database name="db" version="1">
      <genes>
        <gene id="1" protId="tr|P1|P1_HUMAN"/>
        <gene id="2" protId="P2"/>
        <gene id="3" protId="P3"/>
      </genes>
    </database>
  </species>
  <species name="Mus_musculus" NCBITaxId="10090">
    <database name="db" version="1">
      <genes>
        <gene id="4" protId="M1"/>
        <gene id="5" protId="M2"/>
      </genes>
    </database>
  </species>
  <groups>
    <orthologGroup id="HOG:A">
      <paralogGroup>
        <geneRef id="1"/>
        <geneRef id="2"/>
      </paralogGroup>
      <orthologGroup id="HOG:A.1">
        <geneRef id="4"/>
      </orthologGroup>
    </orthologGroup>
    <orthologGroup id="HOG:B">
      <geneRef id="3"/>
      <geneRef id="5"/>
    </orthologGroup>
  </groups>
</orthoXML>
'''


def parsed():
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "x.orthoxml")
        with open(path, "w") as f:
            f.write(XML)
        return ox.parse_orthoxml(path)


def test_every_gene_of_a_top_level_group_is_kept():
    groups, _ = parsed()
    assert list(groups) == ["OG_0000001", "OG_0000002"]
    # paralogs 1 and 2 and the nested ortholog 4 are all in the first group
    assert sorted(groups["OG_0000001"]) == ["1", "2", "4"]
    assert sorted(groups["OG_0000002"]) == ["3", "5"]


def test_species_comes_from_the_species_block():
    _, genes = parsed()
    assert genes["1"]["species"] == "Homo_sapiens"
    assert genes["5"]["species"] == "Mus_musculus"
    assert genes["1"]["protId"] == "tr|P1|P1_HUMAN"
    assert genes["2"]["geneId"] is None


def test_every_gene_in_the_species_blocks_is_recorded():
    _, genes = parsed()
    assert sorted(genes) == ["1", "2", "3", "4", "5"]


def test_namespace_is_stripped_from_tags():
    assert ox._local_tag("{http://orthoXML.org/2011/}species") == "species"
    assert ox._local_tag("species") == "species"


def test_gene_ids_are_made_safe_for_alignment_and_tree_tools():
    assert ox.sanitise_gene_id("tr|B8BSL4|B8BSL4_THAPS") == "tr_B8BSL4_B8BSL4_THAPS"
    assert ox.sanitise_gene_id("a b/c") == "a_b_c"
    assert ox.sanitise_gene_id("keep-this_one1") == "keep-this_one1"
    assert ox.sanitise_gene_id("") == "GENE_UNKNOWN"


def test_orthogroup_names_are_made_safe_for_file_names():
    assert ox.sanitise_og_name("HOG:0001/a") == "HOG_0001_a"
    assert ox.sanitise_og_name("") == "OG_UNKNOWN"


def test_proteomes_are_indexed_per_species_in_sorted_order():
    with tempfile.TemporaryDirectory() as tmp:
        for name, seqs in (("Mus_musculus.fa", ">M1\nMKV\n>M2\nLLA\n"),
                           ("Homo_sapiens.fa", ">P1 desc\nAAA\n"),
                           ("notes.txt", "ignore me")):
            with open(os.path.join(tmp, name), "w") as f:
                f.write(seqs)
        species_to_id, index = ox.index_proteomes(tmp)
    assert species_to_id == {"Homo_sapiens": 1, "Mus_musculus": 2}
    assert sorted(index["Mus_musculus"]) == ["M1", "M2"]
    assert str(index["Homo_sapiens"]["P1"].seq) == "AAA"


def test_the_same_accession_in_two_species_is_allowed():
    with tempfile.TemporaryDirectory() as tmp:
        for name in ("A.fa", "B.fa"):
            with open(os.path.join(tmp, name), "w") as f:
                f.write(">same\nMKV\n")
        _, index = ox.index_proteomes(tmp)
    assert "same" in index["A"] and "same" in index["B"]


def test_a_repeated_accession_inside_one_proteome_is_an_error():
    with tempfile.TemporaryDirectory() as tmp:
        with open(os.path.join(tmp, "A.fa"), "w") as f:
            f.write(">dup\nMKV\n>dup\nLLA\n")
        try:
            ox.index_proteomes(tmp)
        except ValueError as err:
            assert "dup" in str(err)
            return
    raise AssertionError("a repeated accession was accepted")


def test_an_empty_proteome_folder_explains_what_to_do():
    with tempfile.TemporaryDirectory() as tmp:
        os.makedirs(os.path.join(tmp, "proteome"))
        try:
            ox.index_proteomes(tmp)
        except ValueError as err:
            assert "subdirector" in str(err)
            return
    raise AssertionError("an empty proteome folder was accepted")
