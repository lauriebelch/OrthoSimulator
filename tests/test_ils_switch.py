'''
Unit tests for the optional ILS (incomplete lineage sorting) switch.

ILS is off by default. It is switched on only when --ils-theta is set.
When it is off, the coalescent rearrangement step must be skipped and branch
relaxation must read the original pruned tree from SaGePhy.

These tests replace the slow simulation steps with recorders, so they run
without Java, IQ-TREE or any input data.
'''

import os
import sys

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)

import subprocess
import tempfile
from types import SimpleNamespace

import scripts.OrthogroupSimulation.orthogroup_simulation_utils as og


##########################################
######### helpers ########################
##########################################

# the steps that RunOrthogroup calls, in order, before the alignment checks
STEPS_BEFORE_ALIGNMENT = [
    "SimulateGeneTree",
    "RelaxBranchLengths",
    "RelabelRelaxedTree",
    "SelectSequence",
    "MakePartitionPFAM",
    "SimulateAlignmentPartition",
]


def run_orthogroup_and_record_calls(args):
    """
    Run RunOrthogroup with every step replaced by a recorder.
    Returns the names of the steps that were called, in order.
    """
    names = STEPS_BEFORE_ALIGNMENT + [
        "ApplyILSCoalescent",
        "GapShuffleAlignment",
        "CheckAlignmentHealth",
        "CountDuplicationLossTransfer",
        "Orthologs",
        "OrthologParsing",
    ]
    calls = []
    originals = {}

    def make_recorder(name):
        def recorder(*a, **k):
            calls.append(name)
            if name == "SimulateAlignmentPartition":
                return True      # pretend AliSim succeeded
            if name == "CheckAlignmentHealth":
                return False     # pretend the alignment is healthy
            return None
        return recorder

    for name in names:
        originals[name] = getattr(og, name)
        setattr(og, name, make_recorder(name))
    try:
        og.RunOrthogroup("1", args)
    finally:
        for name, func in originals.items():
            setattr(og, name, func)
    return calls


def relax_branch_lengths_and_get_input_tree(args):
    """
    Run RelaxBranchLengths with subprocess.run replaced by a recorder.
    Returns the tree path that was passed to SaGePhy's BranchRelaxer.
    """
    captured = []
    original_run = og.subprocess.run

    def fake_run(cmd, *a, **k):
        captured.append(cmd)
        # BranchRelaxer writes the tree named after -o; pretend it did
        with open(cmd[cmd.index("-o") + 1], "w") as f:
            f.write("();")
        return subprocess.CompletedProcess(cmd, 0, stdout="", stderr="")

    og.subprocess.run = fake_run
    try:
        og.RelaxBranchLengths("1", args)
    finally:
        og.subprocess.run = original_run
    assert len(captured) == 1
    cmd = captured[0]
    # the command is: java -jar <jar> BranchRelaxer <input tree> <model> ...
    return cmd[cmd.index("BranchRelaxer") + 1]


def make_relax_args(tmpdir, **extra):
    os.makedirs(os.path.join(tmpdir, "temporary_files"), exist_ok=True)
    return SimpleNamespace(
        o=tmpdir,
        sagephy_path="sagephy.jar",
        relax_model="ACRY07",
        max_start_rate=1,
        sigma_log_mean=0.0,
        sigma_log_sd=0.1,
        **extra,
    )


##########################################
######### the test functions #############
#########################################

def test_ils_is_off_when_theta_is_not_set():
    assert og.ILSEnabled(SimpleNamespace()) is False
    assert og.ILSEnabled(SimpleNamespace(ils_theta=None)) is False
    assert og.ILSEnabled(SimpleNamespace(ils_theta="")) is False


def test_ils_is_on_when_theta_is_set():
    assert og.ILSEnabled(SimpleNamespace(ils_theta="10")) is True
    assert og.ILSEnabled(SimpleNamespace(ils_theta=10.0)) is True


def test_run_orthogroup_skips_ils_by_default():
    calls = run_orthogroup_and_record_calls(SimpleNamespace(o="unused"))
    assert "ApplyILSCoalescent" not in calls
    # the steps before the alignment checks run in the same order as before ILS existed
    assert calls[:len(STEPS_BEFORE_ALIGNMENT)] == STEPS_BEFORE_ALIGNMENT


def test_run_orthogroup_applies_ils_when_theta_is_set():
    calls = run_orthogroup_and_record_calls(SimpleNamespace(o="unused", ils_theta="10"))
    assert "ApplyILSCoalescent" in calls
    # ILS comes straight after the gene tree is simulated and before branch relaxation
    assert calls.index("SimulateGeneTree") + 1 == calls.index("ApplyILSCoalescent")
    assert calls.index("ApplyILSCoalescent") + 1 == calls.index("RelaxBranchLengths")


def test_branch_relaxation_reads_original_tree_when_ils_is_off():
    with tempfile.TemporaryDirectory() as tmpdir:
        tree = relax_branch_lengths_and_get_input_tree(make_relax_args(tmpdir))
        assert tree.endswith("1.pruned.tree")


def test_branch_relaxation_reads_ils_tree_when_ils_is_on():
    with tempfile.TemporaryDirectory() as tmpdir:
        tree = relax_branch_lengths_and_get_input_tree(
            make_relax_args(tmpdir, ils_theta="10")
        )
        assert tree.endswith("1.pruned_ils.tree")
