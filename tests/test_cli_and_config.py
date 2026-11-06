'''
Unit tests for the command line and config file handling in OrthoSimulator.py.

Parameters can come from a config file, the command line, or both. The
command line wins. Required parameters must be present, unknown config keys
must give a clear message, and an output folder that already exists is
refused so a finished run is never overwritten.
'''

import importlib.util
import os
import sys
import tempfile

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)

spec = importlib.util.spec_from_file_location("OrthoSimulator", os.path.join(REPO_ROOT, "OrthoSimulator.py"))
cli = importlib.util.module_from_spec(spec)
spec.loader.exec_module(cli)

EXAMPLE_CONFIG = os.path.join(REPO_ROOT, "ExampleData", "example_simulation_config.txt")
GA_CONFIG = os.path.join(REPO_ROOT, "ExampleData", "GA_CONFIG.txt")


def all_dests():
    a, b, c, d = cli.args_dicts()
    flags = [v[0] for dct in (a, b, c, d) for v in dct.values()]
    return [f.lstrip("-").replace("-", "_") for f in flags]


def blank_args(**given):
    '''what the argument parser returns: every option present, None if not given'''
    args = {dest: None for dest in all_dests()}
    args.update({"og_sim": False, "GA": False, "get_data": False, "PFAM": False, "threads": 1})
    args.update(given)
    return args


def write(path, text):
    with open(path, "w") as f:
        f.write(text)


def exits(fn):
    try:
        fn()
    except SystemExit:
        return True
    return False


##########################################
######### the option tables ##############
##########################################

def test_every_command_line_flag_is_unique():
    flags = all_dests()
    assert len(flags) == len(set(flags))


def test_every_key_in_the_example_configs_is_a_known_option():
    known = set(all_dests())
    for path in (EXAMPLE_CONFIG, GA_CONFIG):
        with open(path) as f:
            for line in f:
                if line.strip() and not line.startswith("#"):
                    assert line.split("=", 1)[0].strip() in known, (path, line)


def test_the_example_simulation_config_supplies_every_required_simulation_option():
    required = cli.requirements_dicts()["og_sim"]
    with open(EXAMPLE_CONFIG) as f:
        keys = {line.split("=", 1)[0].strip() for line in f if line.strip()}
    # output is a command-line option; everything else must be in the config
    assert set(required) - keys == {"output"}


def test_ils_theta_is_optional_and_not_required():
    assert "ils_theta" not in cli.requirements_dicts()["og_sim"]
    assert "ils_theta" in all_dests()


##########################################
######### check_options ##################
##########################################

def test_missing_required_flags_stop_the_run():
    args = blank_args(og_sim=True, output="x")
    assert exits(lambda: cli.check_options("og_sim", args))


def test_a_config_file_supplies_the_missing_options():
    args = blank_args(og_sim=True, output="x", Config=EXAMPLE_CONFIG)
    params = cli.check_options("og_sim", args)
    assert params["species"] == "Mnemiopsis_leidyi"
    assert params["Orthogroups"] == "100"


def test_the_command_line_wins_over_the_config_file():
    args = blank_args(og_sim=True, output="x", Config=EXAMPLE_CONFIG, Orthogroups="7")
    params = cli.check_options("og_sim", args)
    assert params["Orthogroups"] == "7"


def test_ils_theta_from_the_config_file_reaches_the_run():
    with tempfile.TemporaryDirectory() as tmp:
        cfg = os.path.join(tmp, "c.txt")
        with open(EXAMPLE_CONFIG) as f:
            write(cfg, f.read().rstrip("\n") + "\nils_theta=10\n")
        params = cli.check_options("og_sim", blank_args(og_sim=True, output="x", Config=cfg))
    assert params["ils_theta"] == "10"


def test_ils_is_absent_when_not_given():
    params = cli.check_options("og_sim", blank_args(og_sim=True, output="x", Config=EXAMPLE_CONFIG))
    assert "ils_theta" not in params


def test_a_missing_config_file_stops_the_run():
    args = blank_args(og_sim=True, output="x", Config="/no/such/config.txt")
    assert exits(lambda: cli.check_options("og_sim", args))


def test_a_missing_species_tree_stops_the_run():
    args = blank_args(og_sim=True, output="x", Config=EXAMPLE_CONFIG, ultrametric_tree="/no/such/tree.nwk")
    assert exits(lambda: cli.check_options("og_sim", args))


def test_blank_lines_and_comments_in_a_config_file_are_ignored():
    with tempfile.TemporaryDirectory() as tmp:
        cfg = os.path.join(tmp, "c.txt")
        with open(EXAMPLE_CONFIG) as f:
            write(cfg, "# my settings\n\n" + f.read().rstrip("\n") + "\n\n# end\n")
        params = cli.check_options("og_sim", blank_args(og_sim=True, output="x", Config=cfg))
    assert params["species"] == "Mnemiopsis_leidyi"


def test_a_value_containing_an_equals_sign_is_read_whole():
    with tempfile.TemporaryDirectory() as tmp:
        cfg = os.path.join(tmp, "c.txt")
        with open(EXAMPLE_CONFIG) as f:
            write(cfg, f.read().rstrip("\n") + "\nrelax_model=ACRY07=x\n")
        # relax_model is already in the example config, so the first value wins;
        # use a fresh run without it to see the value passed through whole
        text = open(cfg).read().replace("relax_model=ACRY07\n", "")
        write(cfg, text)
        params = cli.check_options("og_sim", blank_args(og_sim=True, output="x", Config=cfg))
    assert params["relax_model"] == "ACRY07=x"


def test_a_config_line_without_equals_gives_a_clear_message():
    import io
    from contextlib import redirect_stdout
    with tempfile.TemporaryDirectory() as tmp:
        cfg = os.path.join(tmp, "c.txt")
        write(cfg, "species Homo_sapiens\n")
        out = io.StringIO()
        with redirect_stdout(out):
            stopped = exits(lambda: cli.check_options("og_sim", blank_args(og_sim=True, output="x", Config=cfg)))
    assert stopped and "line 1" in out.getvalue()


def test_an_unknown_config_key_gives_a_clear_message():
    import io
    from contextlib import redirect_stdout
    with tempfile.TemporaryDirectory() as tmp:
        cfg = os.path.join(tmp, "c.txt")
        with open(EXAMPLE_CONFIG) as f:
            write(cfg, f.read().rstrip("\n") + "\nnot_an_option=1\n")
        out = io.StringIO()
        with redirect_stdout(out):
            stopped = exits(lambda: cli.check_options("og_sim", blank_args(og_sim=True, output="x", Config=cfg)))
    assert stopped
    assert "not_an_option" in out.getvalue()


##########################################
######### choosing a workflow ############
##########################################

def test_no_workflow_flag_stops_the_run():
    assert exits(lambda: cli.Process_Args(blank_args(), REPO_ROOT))


def test_two_workflow_flags_stop_the_run():
    assert exits(lambda: cli.Process_Args(blank_args(og_sim=True, GA=True), REPO_ROOT))


##########################################
######### output folder ##################
##########################################

def test_an_existing_output_folder_is_refused():
    with tempfile.TemporaryDirectory() as tmp:
        assert exits(lambda: cli.create_outputs_and_configs({"output": tmp}, REPO_ROOT, 1))


def test_the_run_config_records_the_settings_but_not_output_or_config_paths():
    with tempfile.TemporaryDirectory() as tmp:
        out = os.path.join(tmp, "run")
        path, cfg = cli.create_outputs_and_configs(
            {"output": out, "Config": "c.txt", "Orthogroups": "5", "species": "Homo_sapiens"}, REPO_ROOT, 4)
        with open(cfg) as f:
            text = f.read()
    assert "Orthogroups=5" in text and "species=Homo_sapiens" in text
    assert "threads=4" in text
    assert "output=" not in text and "Config=" not in text


def test_the_output_folder_can_be_inside_folders_that_do_not_exist_yet():
    with tempfile.TemporaryDirectory() as tmp:
        out = os.path.join(tmp, "new_parent", "run")
        path, cfg = cli.create_outputs_and_configs({"output": out, "Orthogroups": "5"}, REPO_ROOT, 1)
        assert os.path.isdir(out) and os.path.isfile(cfg)
