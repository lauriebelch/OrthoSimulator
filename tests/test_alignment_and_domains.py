'''
Unit tests for the helpers that read and reshape sequences and alignments:
alignment reading and writing, gap blocks, the gap-position profile,
sequence shuffling, PFAM domain extraction and domain-model selection.
'''

import os
import random
import sys
import tempfile
from types import SimpleNamespace

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)

import scripts.OrthogroupSimulation.orthogroup_simulation_utils as og


def write(path, text):
    with open(path, "w") as f:
        f.write(text)


##########################################
######### alignments #####################
##########################################

def test_alignment_write_then_read_gives_the_same_alignment():
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "a.fa")
        names = ["sp1_1_1", "sp2_1_1", "sp3_1_1"]
        seqs = ["MKV-LL", "MKVALL", "M--ALL"]
        og.WriteAlignment(names, seqs, path)
        assert og.ReadAlignment(path) == (names, seqs)


def test_read_alignment_joins_wrapped_lines_and_skips_blank_lines():
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "a.fa")
        write(path, ">x\nMKV\nLLA\n\n>y\nMK-\nLLA\n")
        assert og.ReadAlignment(path) == (["x", "y"], ["MKVLLA", "MK-LLA"])


def test_column_gap_fractions():
    seqs = ["A-C-", "A-CC", "AACC", "A-C-"]
    # column 1 has gaps in 3 of 4 sequences, column 3 in 2 of 4
    assert og.GetColumnGapFractions(seqs) == [0.0, 0.75, 0.0, 0.5]
    assert og.GetColumnGapFractions([]) == []


def test_gap_blocks_are_runs_of_columns_with_any_gap():
    seqs = ["AA--CC-A",
            "AAAACCCA",
            "AAAACC-A"]
    blocks, fractions = og.GetGapColumnBlocks(seqs)
    assert blocks == [(2, 3), (6, 6)]
    assert len(fractions) == 8


def test_gap_blocks_include_a_block_that_runs_to_the_end():
    blocks, _ = og.GetGapColumnBlocks(["AAA--", "AAAAA"])
    assert blocks == [(3, 4)]


def test_no_gaps_means_no_blocks():
    blocks, _ = og.GetGapColumnBlocks(["AAAA", "AAAA"])
    assert blocks == []


##########################################
######### gap-position profile ###########
##########################################

def profile_text(rows, n_bins=100, total_col="total_gaps"):
    header = ["alignment", total_col] + [f"bin_{i}" for i in range(n_bins)]
    lines = ["\t".join(header)]
    for name, total, probs in rows:
        lines.append("\t".join([name, str(total)] + [str(p) for p in probs]))
    return "\n".join(lines) + "\n"


def test_gap_profile_ranks_bins_from_least_to_most_gap_enriched():
    probs = [0.0] * 100
    probs[10], probs[40], probs[70] = 0.2, 0.5, 0.3
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "p.tsv")
        write(path, profile_text([("aln1", 12, probs)]))
        rows = og.ReadGapPositionProfile(path)
    assert len(rows) == 100
    values = [p for _, p in rows]
    assert values == sorted(values)
    assert [b for b, _ in rows][-3:] == [10, 70, 40]


def test_gap_profile_ignores_alignments_without_gaps():
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "p.tsv")
        write(path, profile_text([("aln1", 0, [0.01] * 100)]))
        assert og.ReadGapPositionProfile(path) == []


def test_gap_profile_rejects_a_table_without_total_gaps():
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "p.tsv")
        write(path, profile_text([("aln1", 5, [0.01] * 100)], total_col="n"))
        try:
            og.ReadGapPositionProfile(path)
        except ValueError as err:
            assert "total_gaps" in str(err)
        else:
            raise AssertionError("a profile without total_gaps was accepted")


def test_gap_profile_rejects_the_wrong_number_of_bins():
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "p.tsv")
        write(path, profile_text([("aln1", 5, [0.01] * 50)], n_bins=50))
        try:
            og.ReadGapPositionProfile(path)
        except ValueError as err:
            assert "100" in str(err)
        else:
            raise AssertionError("a profile with 50 bins was accepted")


def test_the_shipped_example_gap_profile_can_be_read():
    random.seed(0)
    rows = og.ReadGapPositionProfile(os.path.join(REPO_ROOT, "ExampleData", "gap_position_profile_counts.tsv"))
    assert len(rows) == 100


##########################################
######### shuffling ######################
##########################################

def test_shuffle_keeps_length_and_residue_counts():
    random.seed(1)
    seq = "".join(random.choice("ACDEFGHIKLMNPQRSTVWY") for _ in range(237))
    shuffled = og.ShuffleSequence(seq, 20)
    assert len(shuffled) == len(seq)
    assert sorted(shuffled) == sorted(seq)


def test_shuffle_moves_whole_blocks_not_single_residues():
    random.seed(2)
    seq = "".join(chr(65 + i % 20) * 20 for i in range(10))   # ten distinct blocks of 20
    shuffled = og.ShuffleSequence(seq, 20)
    blocks = [shuffled[i:i + 20] for i in range(0, len(shuffled), 20)]
    assert all(len(set(b)) == 1 for b in blocks)


def test_shuffle_of_a_short_sequence_is_unchanged():
    assert og.ShuffleSequence("MKVL", 20) == "MKVL"


##########################################
######### FASTA dictionary ###############
##########################################

def test_fasta_dict_uses_the_first_word_of_the_header_and_joins_lines():
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "g.fa")
        write(path, ">gene1 some description\nMKV\nLLA\n>gene2\nAAA\n")
        d = og.MakeFastaDict(SimpleNamespace(f=path))
    assert d == {"gene1": "MKVLLA", "gene2": "AAA"}


##########################################
######### PFAM domains ###################
##########################################

def scan_line(gene, start, end, hit_type):
    # pfam_scan columns: seq id, aln start, aln end, env start, env end,
    # hmm acc, hmm name, type, ...
    return f"{gene} {start} {end} {start} {end} PF00001.1 Name {hit_type} 1 100 100 50.0 1e-10 1 No_clan\n"


def domains_for(gene, scan_text):
    with tempfile.TemporaryDirectory() as tmp:
        os.makedirs(os.path.join(tmp, "temporary_files"))
        write(os.path.join(tmp, "temporary_files", "1.txt"), f"starting_gene_id={gene}\n")
        scan = os.path.join(tmp, "scan.txt")
        write(scan, scan_text)
        args = SimpleNamespace(o=tmp, PFAM=scan)
        return og.ExtractPfamDomains("1", args)


def test_domains_are_found_only_for_the_starting_gene():
    text = ("# comment\n" + scan_line("geneA", 10, 50, "Domain") + scan_line("geneB", 5, 30, "Domain"))
    hits = domains_for("geneA", text)
    assert [(h["start"], h["end"]) for h in hits] == [(10, 50)]


def test_only_domains_and_motifs_are_used():
    text = (scan_line("geneA", 10, 50, "Domain") + scan_line("geneA", 60, 80, "Family")
            + scan_line("geneA", 90, 99, "Motif") + scan_line("geneA", 100, 120, "Repeat")
            + scan_line("geneA", 130, 140, "Disordered"))
    hits = domains_for("geneA", text)
    assert [(h["start"], h["end"]) for h in hits] == [(10, 50), (90, 99)]


def test_overlapping_domains_are_merged_and_the_result_is_sorted():
    text = (scan_line("geneA", 100, 150, "Domain") + scan_line("geneA", 10, 60, "Domain")
            + scan_line("geneA", 50, 80, "Domain"))
    hits = domains_for("geneA", text)
    assert [(h["start"], h["end"]) for h in hits] == [(10, 80), (100, 150)]


def test_a_gene_with_no_domains_gives_an_empty_list():
    assert domains_for("geneA", scan_line("geneB", 1, 10, "Domain")) == []


def test_random_domain_model_is_taken_from_the_alisim_model_column():
    random.seed(0)
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "m.csv")
        write(path, "domain,alisim_model\nPF1,LG+G\nPF2,\nPF3,JTT+G\n")
        picks = {og.SelectRandomDomainModel(SimpleNamespace(d=path)) for _ in range(50)}
    assert picks == {"LG+G", "JTT+G"}      # the row with no model is never picked
