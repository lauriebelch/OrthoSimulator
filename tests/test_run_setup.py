'''
Unit tests for the setup code that runs before a simulation starts
(scripts/OrthogroupSimulation/simulate_orthogroup_main.py): config reading,
path handling, and the checks that stop a bad run before any work is done.
'''

import os
import sys
import tempfile
from types import SimpleNamespace

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)

import scripts.OrthogroupSimulation.orthogroup_simulation_utils as og
import scripts.OrthogroupSimulation.simulate_orthogroup_main as sm


def write(path, text):
    with open(path, "w") as f:
        f.write(text)


def raises(exc, fn):
    try:
        fn()
    except exc as err:
        return str(err)
    raise AssertionError(f"{exc.__name__} was not raised")


##########################################
######### config helpers #################
##########################################

def test_read_config_skips_blank_and_comment_lines():
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "c.txt")
        write(path, "# comment\n\na=1\nb = two \n")
        assert sm.read_config(path) == {"a": "1", "b": "two"}


def test_read_config_keeps_equals_signs_inside_values():
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "c.txt")
        write(path, "path=/a/b=c/d\n")
        assert sm.read_config(path) == {"path": "/a/b=c/d"}


def test_read_config_rejects_a_line_without_equals():
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "c.txt")
        write(path, "oops\n")
        assert "line 1" in raises(ValueError, lambda: sm.read_config(path))


def test_require_config_returns_the_value_or_names_the_missing_key():
    assert sm.require_config({"a": "1"}, "a") == "1"
    assert "species" in raises(ValueError, lambda: sm.require_config({"a": "1"}, "species"))
    assert "species" in raises(ValueError, lambda: sm.require_config({"species": "  "}, "species"))


def test_relative_paths_are_resolved_from_the_project_root():
    assert sm.resolve_path("ExampleData/x.nwk", "/proj") == "/proj/ExampleData/x.nwk"
    assert sm.resolve_path("/abs/x.nwk", "/proj") == "/abs/x.nwk"


def test_first_existing_file_returns_the_first_match_and_lists_the_tries():
    with tempfile.TemporaryDirectory() as tmp:
        real = os.path.join(tmp, "real.fa")
        write(real, "x")
        assert sm.first_existing_file([None, os.path.join(tmp, "no.fa"), real], "genome") == real
        message = raises(ValueError, lambda: sm.first_existing_file([os.path.join(tmp, "no.fa")], "genome"))
    assert "genome" in message and "no.fa" in message


def test_the_project_root_is_the_repository_folder():
    assert os.path.samefile(sm.infer_project_root({}), REPO_ROOT)


##########################################
######### parameter files ################
##########################################

def test_parameter_file_values_get_the_right_types():
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "p.txt")
        write(path, "n=100\nrate=0.5\nsmall=1e-3\nname=ACRY07\nflag=true\nother=off\n\nbad line\n")
        p = og.LoadParameters(path)
    assert p["n"] == 100 and isinstance(p["n"], int)
    assert p["rate"] == 0.5 and p["small"] == 1e-3
    assert p["name"] == "ACRY07"
    assert p["flag"] is True and p["other"] is False
    assert "bad line" not in p


def test_a_negative_number_is_read_as_a_number():
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "p.txt")
        write(path, "x=-1.2226\ny=-3\n")
        p = og.LoadParameters(path)
    assert p["x"] == -1.2226 and p["y"] == -3


##########################################
######### the ILS check ##################
##########################################

def make_args(tmp, **extra):
    paths = {}
    for name in ("config", "tree", "gap", "param", "genome", "scan", "csv", "jar", "iqtree"):
        paths[name] = os.path.join(tmp, name)
        write(paths[name], "x")
    write(paths["tree"], "((A:1,B:1):1,C:2);")
    os.chmod(paths["iqtree"], 0o755)
    write(paths["param"], "".join(f"{k}={v}\n" for k, v in {
        "sagephy_path": paths["jar"], "iqtree_path": paths["iqtree"], **extra}.items()))
    return SimpleNamespace(config=paths["config"], s=paths["tree"], gap_profile=paths["gap"],
                           p=paths["param"], f=paths["genome"], PFAM=paths["scan"], d=paths["csv"])


def test_a_good_run_passes_the_input_check_without_ils():
    with tempfile.TemporaryDirectory() as tmp:
        args = make_args(tmp)
        sm.check_inputs_and_load_parameters(args)
        assert not og.ILSEnabled(args)


def test_zero_and_negative_theta_are_refused_with_a_hint():
    for bad in ("0", "-1", "-0.5"):
        with tempfile.TemporaryDirectory() as tmp:
            args = make_args(tmp, ils_theta=bad)
            message = raises(ValueError, lambda: sm.check_inputs_and_load_parameters(args))
            assert "ils_theta" in message and "leave out" in message


def test_positive_theta_is_accepted_and_switches_ils_on():
    with tempfile.TemporaryDirectory() as tmp:
        args = make_args(tmp, ils_theta="10")
        sm.check_inputs_and_load_parameters(args)
        assert og.ILSEnabled(args)


def test_a_missing_input_file_is_named_in_the_error():
    with tempfile.TemporaryDirectory() as tmp:
        args = make_args(tmp)
        os.remove(args.s)
        assert "species tree" in raises(ValueError, lambda: sm.check_inputs_and_load_parameters(args))


def test_a_non_executable_iqtree_is_reported():
    with tempfile.TemporaryDirectory() as tmp:
        args = make_args(tmp)
        os.chmod(os.path.join(tmp, "iqtree"), 0o644)
        assert "not executable" in raises(ValueError, lambda: sm.check_inputs_and_load_parameters(args))


def test_zero_or_negative_gbc_is_refused_before_any_simulation_starts():
    for bad in ("0", "-2"):
        with tempfile.TemporaryDirectory() as tmp:
            args = make_args(tmp, gbc=bad)
            assert "gbc must be > 0" in raises(ValueError, lambda: sm.check_inputs_and_load_parameters(args))


def test_a_positive_gbc_is_accepted():
    with tempfile.TemporaryDirectory() as tmp:
        sm.check_inputs_and_load_parameters(make_args(tmp, gbc="6.75"))


def test_a_missing_sagephy_tree_is_reported_with_sagephys_message():
    import subprocess
    original = og.subprocess.run

    def fake_run(cmd, *a, **k):
        return subprocess.CompletedProcess(cmd, 0, stdout="Unknown option -gbc 0\n", stderr="")

    og.subprocess.run = fake_run
    try:
        with tempfile.TemporaryDirectory() as tmp:
            os.makedirs(os.path.join(tmp, "temporary_files"))
            args = SimpleNamespace(o=tmp, sagephy_path="x.jar", gbc=0, leaf_sampling_probability=1,
                                   replacement_prob=0, s="tree.nwk", max_duplication_rate=0.5,
                                   max_loss_rate=0.5, max_transfer_rate=0)
            message = raises(RuntimeError, lambda: og.SimulateGeneTree("3", args))
    finally:
        og.subprocess.run = original
    assert "orthogroup 3" in message and "Unknown option" in message


def test_an_ultrametric_species_tree_is_accepted_and_rounding_is_tolerated():
    with tempfile.TemporaryDirectory() as tmp:
        for newick in ("((A:1,B:1):1,C:2);", "((A:1.000001,B:1):1,C:2.0000005);"):
            path = os.path.join(tmp, "t.nwk")
            write(path, newick)
            sm.check_species_tree_is_ultrametric(path)


def test_a_species_tree_that_is_not_ultrametric_is_refused():
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "t.nwk")
        write(path, "((A:1,B:1):1,C:3);")
        assert "not ultrametric" in raises(ValueError, lambda: sm.check_species_tree_is_ultrametric(path))


def test_a_species_tree_without_branch_lengths_is_refused():
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "t.nwk")
        write(path, "((A,B),C);")
        assert "branch lengths" in raises(ValueError, lambda: sm.check_species_tree_is_ultrametric(path))


def test_the_shipped_example_species_tree_is_ultrametric():
    sm.check_species_tree_is_ultrametric(os.path.join(REPO_ROOT, "ExampleData", "species_tree.nwk"))


def relax_with_fake_relaxer(tmp, succeed_on_attempt, **extra):
    '''run RelaxBranchLengths with a fake BranchRelaxer that writes its tree only
    from the given attempt onwards. Returns the sigma values it was called with.'''
    import subprocess
    sigmas = []

    def fake_run(cmd, *a, **k):
        sigmas.append(float(cmd[cmd.index("BranchRelaxer") + 4]))
        if len(sigmas) >= succeed_on_attempt:
            with open(cmd[cmd.index("-o") + 1], "w") as f:
                f.write("();")
        return subprocess.CompletedProcess(cmd, 0, stdout="Failed to create valid rates", stderr="")

    original = og.subprocess.run
    og.subprocess.run = fake_run
    try:
        os.makedirs(os.path.join(tmp, "temporary_files"), exist_ok=True)
        args = SimpleNamespace(o=tmp, sagephy_path="x.jar", relax_model="ACRY07", max_start_rate=1,
                               sigma_log_mean=0.0, sigma_log_sd=0.5, **extra)
        og.RelaxBranchLengths("1", args)
    finally:
        og.subprocess.run = original
    return sigmas


def test_branch_relaxation_tries_again_with_a_new_sigma_when_no_tree_is_written():
    with tempfile.TemporaryDirectory() as tmp:
        sigmas = relax_with_fake_relaxer(tmp, succeed_on_attempt=3)
    assert len(sigmas) == 3 and len(set(sigmas)) == 3


def test_branch_relaxation_gives_up_with_a_clear_message_after_ten_attempts():
    with tempfile.TemporaryDirectory() as tmp:
        message = raises(RuntimeError, lambda: relax_with_fake_relaxer(tmp, succeed_on_attempt=999))
        assert "10 attempts" in message and "orthogroup 1" in message and "ils_theta" not in message
        message = raises(RuntimeError, lambda: relax_with_fake_relaxer(tmp, succeed_on_attempt=999, ils_theta="0.01"))
        assert "ils_theta" in message


def test_the_sigma_is_capped_at_100():
    with tempfile.TemporaryDirectory() as tmp:
        import subprocess
        sigmas = []
        original = og.subprocess.run

        def fake_run(cmd, *a, **k):
            sigmas.append(float(cmd[cmd.index("BranchRelaxer") + 4]))
            with open(cmd[cmd.index("-o") + 1], "w") as f:
                f.write("();")
            return subprocess.CompletedProcess(cmd, 0, stdout="", stderr="")

        og.subprocess.run = fake_run
        try:
            os.makedirs(os.path.join(tmp, "temporary_files"))
            args = SimpleNamespace(o=tmp, sagephy_path="x.jar", relax_model="ACRY07", max_start_rate=1,
                                   sigma_log_mean=10.0, sigma_log_sd=0.1)
            og.RelaxBranchLengths("1", args)
        finally:
            og.subprocess.run = original
    assert sigmas == [100.0]
