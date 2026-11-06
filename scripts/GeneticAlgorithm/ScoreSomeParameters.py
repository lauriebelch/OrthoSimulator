 #!/usr/bin/env python3

import argparse
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
import subprocess
import uuid
import shutil
import argparse
#SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))


def run_orthogroup_simulation(
    complete_parameters,
    output_abolsute_path,
    threads,
    current_file_path,
    config_file_path,):

    import os
    import time
    import multiprocessing as mp
    from pathlib import Path
    from types import SimpleNamespace
    from importlib.machinery import SourceFileLoader

    # ------------ setup -----------
    ## threads
    threads = int(threads)
    if threads < 1:
        raise ValueError("threads must be >= 1")

    ## config and pathing
    config_file_path = os.path.abspath(config_file_path)
    output_abolsute_path = os.path.abspath(output_abolsute_path)
    scripts = os.path.join(output_abolsute_path,"sim")	

    import scripts.OrthogroupSimulation.orthogroup_simulation_utils as og
    import scripts.OrthogroupSimulation.simulate_orthogroup_main as og_main


    if not os.path.isfile(config_file_path):
        raise ValueError(f"OrthoSimulator config not found: {config_file_path}")

    if not os.path.isdir(output_abolsute_path):
        raise ValueError(
            f"Output directory not found: {output_abolsute_path}\n"
            "OrthoSimulator.py should create this directory before running the simulation."
        )

    config = og_main.read_config(config_file_path)
    ##args
    args = SimpleNamespace(
        config=config_file_path,
        output=output_abolsute_path,
        threads=threads,
    )
    args = og_main.configure_legacy_args(args, config)

    og_main.create_subfolders(args.o)
    og_main.check_inputs_and_load_parameters(args)
    og.THREADS = threads
    #og_main.print_detected_inputs(args)

    ## multiprocesing
    mp.set_start_method("spawn", force=True)
    og.RunOrthogroupMultiProc(
        n=args.n,
        outdir=args.o,
        threads=threads,
        args=args,
    )

    # printing and saving output
    # ------------ workflow -----------

    og.CopyAlignments(args)
    og.CopyTrees(args)
    og.BuildProteomes(args)
    og.ConcatOrthologs(args)
    og.SaveOrthogroups(args)
    og.BuildOrthogroupStats(args)
    og.WriteRunLog(args)
	



def read_config(config_path):
    """
    Read key=value config file.

    Blank lines and comment lines are ignored.
    """
    config = {}

    with open(config_path) as handle:
        for line_number, line in enumerate(handle, start=1):
            line = line.strip()

            if not line or line.startswith("#"):
                continue

            if "=" not in line:
                raise ValueError(
                    f"Bad config line {line_number} in {config_path}: {line}"
                )

            key, value = line.split("=", 1)
            config[key.strip()] = value.strip()

    return config


def modify_param_file(complete_parameters,param_file,current_file_path,default_species_pick,run_dir, dup_rate, loss_rate, gbc_val):
    new_param_file = os.path.join(run_dir,"parameter_file.txt")
    with open(param_file) as f:
        lines = f.readlines()

    ######### needs...
    """
    ultrametric_tree=Training_FastOMA/TrainingResults/SimulationInputs/species_tree.nwk
    species=Mnemiopsis_leidyi
    Orthogroups=100
    threads=4
    """

    with open(new_param_file, "w") as f:
        tree_path = "ultrametric_tree=" + os.path.join(complete_parameters["Orthogroup_train_results"],"TrainingResults","SimulationInputs","species_tree.nwk")
        species = "species="+ default_species_pick
        Orthogroups  = "Orthogroups=" + complete_parameters['Orthogroups']
        gap_file="gap_file="+os.path.join(complete_parameters["Orthogroup_train_results"],"TrainingResults","SimulationInputs","gap_position_profile_counts.tsv")
        threads = "threads=1"
        output = "SimPath=" + run_dir

        f.write(f"{tree_path}\n")
        f.write(f"{species}\n")
        f.write(f"{Orthogroups}\n")
        f.write(f"{threads}\n")
        f.write(f"{output }\n")
        f.write(f"{gap_file }\n")

        for line in lines:
            if line.startswith("max_duplication_rate="):
                f.write(f"max_duplication_rate={dup_rate}\n")
            elif line.startswith("max_loss_rate="):
                f.write(f"max_loss_rate={loss_rate}\n")
            elif line.startswith("gbc="):
                f.write(f"gbc={gbc_val}\n")
            elif line.startswith("sagephy_path"):
                new_sagephy_path = os.path.join(run_dir,"bin","sagephy-1.0.0.jar")
                new_sagephy_path_command = "sagephy_path=" + new_sagephy_path
                f.write(new_sagephy_path_command+"\n")

            elif line.startswith("iqtree"):
                new_iqtree_path = os.path.join(run_dir,"bin","iqtree3")
                new_iqtree_path_command = "iqtree_path=" + new_iqtree_path
                f.write(new_iqtree_path_command+"\n")

            else:
                f.write(line)



def main(dup_rate,loss_rate,gbc_val,complete_parameters,current_file_path,default_species_pick):
    #parser = argparse.ArgumentParser()
    #parser.add_argument("-d", type=float, required=True)
    #parser.add_argument("-l", type=float, required=True)
    #parser.add_argument("-gbc", type=float, required=True)

    #args = parser.parse_args()

    #dup_rate = args.d
    #loss_rate = args.l
    #gbc_val = args.gbc


    # --- create run folder ---
    base_out = os.path.abspath(os.path.join(complete_parameters["output"],"Simulation_Temp_Files"))
    os.makedirs(base_out, exist_ok=True)

    run_dir = os.path.join(base_out, uuid.uuid4().hex)

    os.makedirs(run_dir)
    # -------------- TEMP ------------------
    ## copy over all default_species_related files to the run_dir
    PFAM_Path = os.path.join(current_file_path,"PFAM")
    Copy_PFAM_Path = os.path.join(run_dir,"PFAM")	

    ## need better way...
    #shutil.copytree(PFAM_Path,
    #            Copy_PFAM_Path)


    # --- modify parameter file INSIDE copied folder ---
    param_file = os.path.join(complete_parameters["Orthogroup_train_results"],"TrainingResults","SimulationInputs", "simulation_parameters.txt")
    modify_param_file(complete_parameters,param_file,current_file_path,default_species_pick,run_dir, dup_rate, loss_rate, gbc_val)


    # --- paths to scripts ---
    #simulate_script = os.path.join(SCRIPT_DIR, "simulate_orthogroup.py")
    #stats_script = os.path.join(SCRIPT_DIR, "GetSimulationStats.py")
    #compare_script = os.path.join(SCRIPT_DIR, "CompareDistributions.py")

    # --- copy binary files to new location ---
    binary_bin = os.path.join(current_file_path,"bin")
    run_bin = os.path.join(run_dir,"bin")
    shutil.copytree(binary_bin,run_bin)

    # --- copy sim scripts ---
    #print(current_file_path)
    sim_bin = os.path.join(current_file_path,"scripts","OrthogroupSimulation")
    run_sim = os.path.join(run_dir,"sim")
    shutil.copytree(sim_bin,run_sim )





    #sim_out = os.path.join(run_dir, "sim")
	
    # --- run simulation ---
    #cmd_sim = [
    #    "python", simulate_script,
    #    "--folder", sim_inputs_copy,
    #    "-n", "500",
    #    "-o", sim_out,
    #    "-d", "Mnemiopsis_leidyi"
    #]
    #subprocess.run(cmd_sim, check=True, cwd=SCRIPT_DIR)

    # --- run simulation ---
    config_file_params = read_config(os.path.join(run_dir, "parameter_file.txt"))
    ###################### just need to fix this part now...	


    run_orthogroup_simulation(
        config_file_params,
        run_dir,
        1,
        current_file_path,
        os.path.join(run_dir,"parameter_file.txt"),
    )

    # --- stats ---    old... #stats_script = os.path.join(SCRIPT_DIR, "GetSimulationStats.py")
	
    #cmd_stats = [
    #    "python", stats_script,
    #    "-f", sim_out
    #]
    #
    #subprocess.run(cmd_stats, check=True, cwd=SCRIPT_DIR)
    import scripts.GeneticAlgorithm.GetSimulationStats as GetSimulationStats
    GetSimulationStats.GetSimulationStats_main(run_dir)


    # --- compare ---
    #cmd_compare = [
    #    "python", compare_script,
    #    "-s", sim_out,
    #    "-e", TRAINING_ANALYSIS
    #]

    import scripts.GeneticAlgorithm.CompareDistributions as CompareDistributions
    emp_data_path = os.path.join(complete_parameters["Orthogroup_train_results"],"TrainingResults","Orthogroup_analysis")	
    scores = CompareDistributions.compare_distributions(emp_data_path,run_dir)
    #print(scores)
    #quit()


    #result = subprocess.run(
    #    cmd_compare,
    #    capture_output=True,
    #    text=True,
    #    check=True,
    #    cwd=SCRIPT_DIR
    #)

    # --- extract score ---
    #score = result.stdout.strip().splitlines()[-1]

    return scores,run_dir


if __name__ == "__main__":
    main()
