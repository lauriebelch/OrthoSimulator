'''
Unit tests for OrthogroupSimulation/orthogroup_simulation_utils.py

protocol is to:
1) import a function from the script
2) define an input and output for a case where we know the answer
3) assert that the function returns that answer

test function names must begin with 'test_'
'''

import os
import sys

# make sure the repo root is importable
REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)

# do general imports
import tempfile
from pathlib import Path
from types import SimpleNamespace

# import functions that we want to test
from scripts.OrthogroupSimulation.orthogroup_simulation_utils import (
    SpeciesFromGene,
    CountDuplicationLossTransfer
)

##########################################
######### the test functions #############
#########################################

# SpeciesFromGene

def test_species_from_gene_basic():
    # this function extracts the species name from a simulated gene name
    gene = "homo_sapiens_12_2"
    assert SpeciesFromGene(gene) == "homo_sapiens"

# CountDuplicationLossTransfer

def test_count_duplication_loss_transfer_writes_expected_log_lines():
    # this function counts duplications, losses etc. from files, to be written to logs and summaries
    # write a temp file with a certain number of duplications losses etc.
    with tempfile.TemporaryDirectory() as tmp:
        args = SimpleNamespace(o=tmp)
        tmp_files = Path(tmp) / "temporary_files"
        tmp_files.mkdir()
        (tmp_files / "1.pruned.info").write_text(
            "No. of duplications: 3\n"
            "No. of additive transfers: 1\n"
            "No. of replacing transfers: 2\n"
        )
        (tmp_files / "1.unpruned.info").write_text("No. of losses: 5\n")
        # do the function that is supposed to count these events
        CountDuplicationLossTransfer("1", args)
        log_content = (tmp_files / "1.txt").read_text()
    # assert the values that we expect
    assert "num_losses=5" in log_content
    assert "num_duplications=3" in log_content
    assert "num_additive_transfers=1" in log_content
    assert "num_replacing_transfers=2" in log_content
    assert "total_transfers=3" in log_content  # additive + replacing

# continue with more test functions
