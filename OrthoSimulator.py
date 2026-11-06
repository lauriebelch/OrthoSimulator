# -*- coding: utf-8 -*-
"""
@author: JH & LB
"""

#Simulations workflow wrapper

import argparse
import contextlib
import json
import os
import shutil
import sys
import time
from pathlib import Path

from collections import Counter

# Must be set before numpy import.
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

def path_test(args):
    test = [i.lstrip('-').replace("-","_") for i in ["--ultrametric-tree","--Proteome","--parameter-folder","--tools-proteomes","--tool-output","--tool-tree"]]
    s = 0
    for key,value in args.items():
        if key in test:
            if os.path.isfile(os.path.abspath(value)) == False and os.path.isdir(os.path.abspath(value)) == False:
                print("path not found %s" % (key + ":" + value))
                s = s + 1
    supported_tools = ["OF3","FASTOMA","BROCCOLI","SONICPARANOID2","ORTHOXML"]
    test = [i.lstrip('-').replace("-","_") for i in ["--tool"]]     
    for key,value in args.items():
        if key in supported_tools:
            value_ = value.lstrip().rstrip()
            if value_ not in supported_tools:
                s = s + 1
                print("%s is not a supported tool" % value)
    if s != 0:
        sys.exit()
            

def args_dicts():
    args = {"output": ["--output","OUTPUT"],
            "Input_Tree": ["--ultrametric-tree","TREE_PATH"],
            "Species": ["--species","SPECIES_NAME"],
            "Proteome": ["--Proteome","PATH_TO_PROTEOME"]}
            
    ##### Orthogroups simulations             
    args_og_sim = {
        "prop_invar_mean":["--prop-invar-mean","float"],
        "prop_invar_sd":["--prop-invar-sd","float"],
        "gamma_shape_mean":["--gamma-shape-mean","float"],
        "gamma_shape_sd":["--gamma-shape-sd","float"],
        "max_indel_insert":["--max-indel-insert","float"],
        "max_indel_delete":["--max-indel-delete","float"],
        "indel_size":["--indel-size","float"],
        "max_duplication_rate":["--max-duplication-rate","float"],
        "max_loss_rate":["--max-loss-rate","float"],
        "max_transfer_rate":["--max-transfer-rate","float"],
        "replacement_prob":["--replacement-prob","float"],
        "leaf_sampling_probability":["--leaf-sampling-probability","float"],
        "relax_model":["--relax-model","str"],
        "max_start_rate":["--max-start-rate","float"],
        "sigma_log_mean":["--sigma-log-mean","float"],
        "sigma_log_sd":["--sigma-log-sd","float"],
        "gbc":["--gbc","int"],
        "gap_file":["--gap-file","PATH_TO_GAP_PROFILE"],
    }
            
    ################################################# GA
    args_GA = {
                "num_generations":["--generations","int"],
                "sample_size":["--sample-size","int"],
                "num_parents_mating":["--num-parents-mating","int"],
                "sol_per_pop":["--sol-per-pop","int"],
                "num_genes":["--num-genes","int"],
                "parent_selection_type":["--parent-selection-type","str"],
                "keep_parents":["--keep-parents","bool"],
                "crossover_type":["--crossover-type","str"],
                "crossover_probability":["--crossover-probability","int"],
                "mutation_type":["--mutation-type","str"],
                "mutation_percent_genes":["--mutation-percent-genes","int"],
                "save_solutions":["--save-solutions","bool"],
                "parallel_processing":["--parallel-processing","int"],
                "Orthogroup_train_results":["--Orthogroup-train-results","Path to an orthogroup train folder see --Get-Parameters"]
        
        }
                 
    ################################################# SHARED
    args_shared = {
            "OG_num":["--Orthogroups","NUMBER_OF_ORTHOGROUPS"],            
            "tool":["--tool","OF3,FASTOMA,BROCCOLI,SONICPARANOID2,ORTHOXML"],
            "Tool-output":["--tool-output","PATH_TO_TOOL_OUTPUT (for ORTHOXML: path to the .orthoxml file itself)"],
            "Tool-tree":["--tool-tree","PATH_TO_TOOL_TREE required for: BROCCOLI, SONICPARANOID2, ORTHOXML (these do not produce/carry a usable species tree natively); optional override for FASTOMA"],
            "Tool_input_Proteomes":["--tools-proteomes","PATH_TO_TOOLS_PROTEOMES required for: FASTOMA, BROCCOLI, SONICPARANOID2, ORTHOXML"],
            #"GA_scores":["--parameter-folder","PATH_TO_PARAMETERS note: output of --Get-Parameters"],
            "Config":["--Config","PATH_TO_CONFIG FILE"],
            # NOTE: deliberately NOT part of args_og_sim - it is optional (ILS is
            # off unless it is set) and so must not be in the "og_sim"
            # required-flags list. It is also not GA-tunable and will never
            # appear in a GA-produced config.txt.
            "ils_theta":["--ils-theta","float, optional. Switches on incomplete lineage sorting (ILS). Tree depth in coalescent units (higher = less ILS), must be > 0. If not set, ILS is off. Not tuned by --GA: set it directly."],
            }

    return args,args_og_sim,args_GA,args_shared

def arg_parse_reformat(l):
    return [arg.lstrip('-').replace("-","_") for arg,h in l]

def requirements_dicts():
    args,args_og_sim,args_GA,args_shared = args_dicts()
    requirement_dict = {
        "og_sim" : arg_parse_reformat([args["output"],args["Input_Tree"],args["Species"],args_shared["OG_num"]] + list(args_og_sim.values())),
        "GA" : arg_parse_reformat([args["output"],args_shared["OG_num"]] +list(args_GA.values())), 
        "get_data" : arg_parse_reformat([args["output"],args_shared["tool"],args_shared["Tool-output"]]),
        "PFAM" : arg_parse_reformat([args["output"],args["Proteome"]])}
    return requirement_dict
                     
def Programme_Call():
    parser = argparse.ArgumentParser(
        description=
"""
Orthogroup simulations Toolkit
-----------------------------
Programmes:
Orthogroup simulation
Parameter Tuning
PFAM run and collection
Provide Data for parameters : thing about it : supported options = OF3, FastOMA. 
""",formatter_class=argparse.RawTextHelpFormatter
    )

    args,args_og_sim,args_GA,args_shared = args_dicts()
    
    og_sim_help_list = [args["output"],args["Input_Tree"],args["Species"],args_shared["Config"],args_shared["OG_num"],args_shared["ils_theta"]] + list(args_og_sim.values())
    og_sim_str = [" ".join(i) for i in og_sim_help_list]
    parser.add_argument(
        "--Orthogroup-simulation",
        dest = "og_sim",
        action='store_true',
        default=False,
        help=(
             "\n".join(og_sim_str) + "\n "
        )
    )

    GA_help_list  = [args["output"]] +list(args_GA.values())
    GA_help = [" ".join(i) for i in GA_help_list]
    parser.add_argument(
        "--GA",
        dest = "GA",
        action='store_true',
        default=False,
        help=(
            "\n" +"\n".join(GA_help) + "\n "
        )
    )
    
    
    get_data_from_run_help = [args["output"],args_shared["tool"],args_shared["Tool-output"],args_shared["Tool-tree"],args_shared["Tool_input_Proteomes"]]
    get_data_from_run_help_str = [" ".join(i) for i in get_data_from_run_help]    
    parser.add_argument(
        "--Get-Parameters",
        dest = "get_data",
        action='store_true',
        default=False,
        help=(
            "\n" +"\n".join(get_data_from_run_help_str) + "\n "
        )
    )
    
    pfam_help = [args["output"],args["Proteome"]]
    pfam_help_str =  [" ".join(i) for i in pfam_help]  
    
    parser.add_argument(
        "--PFAM",
        dest = "PFAM",
        action='store_true',
        default=False,
        help=(
            "\n" +"\n".join(pfam_help_str) + "\n "
        )
    )   

    ###################################################################### REAL FUNCTION CALLS
    
    for arg_dict in [args,args_og_sim,args_GA,args_shared]:
        for key, value in arg_dict.items():   

            if key in ["Config"]: ### optional flags
                name, description = value
                parser.add_argument(
                    name,
                    dest = name.lstrip('-').replace("-","_"),
                    default=None,
                    help=argparse.SUPPRESS
                )
                
            else:
                name, description = value
                parser.add_argument(
                    name,
                    dest = name.lstrip('-').replace("-","_"),
                    default=None,
                    help=argparse.SUPPRESS
                )

     

     ###################################################################### thread call        
    parser.add_argument(
        "--threads",
        dest = "threads",
        default=1,
        help=(
            "Threads to be passed to each stage of the workflow"
        )
    )
    
    return parser.parse_args()

def check_missing(requirements,args):
    # "requirements" only tells us what's required for THIS function call.
    # complete_parameters should carry through every flag the user actually
    # supplied (e.g. optional flags like --tools-proteomes/--tool-tree), not
    # just the required ones, otherwise they get silently dropped here and
    # never reach the workflow functions below.
    args_,args_og_sim,args_GA,args_shared = args_dicts()
    known_keys = set(
        arg_parse_reformat(list(args_.values()))
        + arg_parse_reformat(list(args_og_sim.values()))
        + arg_parse_reformat(list(args_GA.values()))
        + arg_parse_reformat(list(args_shared.values()))
    )
    missing = [r for r in requirements if args[r] == None]
    complete_parameters = {k: v for k, v in args.items() if k in known_keys and v is not None}
    return missing,complete_parameters


def check_options(function_call,args):

    ########### check if config is called...
    #args,args_og_sim,args_GA,args_shared = args_dicts()
    requirements  = requirements_dicts()[function_call]

    ############ if config is not called...
    if args['Config'] == None:
        missing,complete_parameters  = check_missing(requirements,args)
        if len(missing) != 0:
            print("The following flags are missing:\n" + ", ".join(missing))
            sys.exit()
            
    elif args['Config'] != None: 
        ## test if i can find the file...
        if os.path.isfile(os.path.abspath(args['Config'])) == False:
            print("Could not find config file at %s" % args['Config'])
            sys.exit()
            
        else:
            ##open and read config...
            with open(os.path.abspath(args['Config'])) as config:
                for line_number, line in enumerate(config, start=1):
                    line = line.strip()
                    # blank lines and comments are allowed
                    if not line or line.startswith("#"):
                        continue
                    if "=" not in line:
                        print("Config file line %d is not in key=value form: %s" % (line_number, line))
                        sys.exit()
                    arg_,value_ = line.split("=", 1)
                    arg = arg_.strip()
                    value = value_.strip()
                    if arg not in args:
                        print("Unknown option in config file (line %d): %s" % (line_number, arg))
                        sys.exit()
                    if args[arg] == None:
                        args[arg] = value
            missing,complete_parameters   = check_missing(requirements,args)
            if len(missing) != 0:
                print("The following flags are missing:\n" + ", ".join(missing))
                sys.exit()
            #### check missing 
            
    ## tested all required paths..
    path_test(complete_parameters)
    ## return complete parameters
    return complete_parameters
    ###
    
def create_outputs_and_configs(complete_parameters,current_file_path, threads):
    ## if results folder already exisits cancel...
    ## return path to results folder path to config 
    #print(complete_parameters)
    output_abolsute_path = os.path.abspath(complete_parameters['output'])
    if os.path.isdir(output_abolsute_path):
        print("Output Directory Path already exists at:\n%s\nPlease specifiy a new path or remove the results directory" % output_abolsute_path)
        sys.exit()
    else:
        try:
            os.makedirs(output_abolsute_path)
        except:
            print("unknown error while creating %s" % output_abolsute_path)
            sys.exit()
            
    config_file_path = os.path.join(output_abolsute_path,"config.txt")
    with open(config_file_path,"a") as config:
        for parameter, value in complete_parameters.items():
            if parameter not in ["output","Config","threads"]:
                config.write("%s=%s\n" % (parameter,value))
        config.write("threads=%s\n" % threads)
        config.write("SimPath=%s" % current_file_path)
    print("created output folder and config file")
    return output_abolsute_path,config_file_path
    
def Process_Args(args,current_file_path):
    
    ## identify function call
    function_calls = {"og_sim":False,"GA":False,"get_data":False,"PFAM":False}
    function_call = None
    for function in function_calls.keys():
        if args[function] == True:
            function_calls[function] = True
            function_call = function
    ## check correct numbers are called.          
    function_call_counts = dict(Counter(function_calls.values()))

    if function_call_counts[False] == 4:
        print("no option selected")
        print("you can see options by running the -h flag")
        sys.exit()
    if function_call_counts[True] > 1:
        print("to many options")
        sys.exit()    
    threads = args['threads']
    complete_parameters = check_options(function_call,args)
    output_abolsute_path,config_file_path = create_outputs_and_configs(complete_parameters,current_file_path,threads)
    
    ##########################
    ### data to pass back to main commands..
    return function_call,output_abolsute_path,config_file_path,threads,complete_parameters
    
    # going to return function call and path containing config file...
    #os.path.abspath("path")

## to run orthogroup simulation
def run_orthogroup_simulation(
    complete_parameters,
    output_abolsute_path,
    threads,
    current_file_path,
    config_file_path
):
    import scripts.OrthogroupSimulation.simulate_orthogroup_main as SimMain

    SimMain.run_orthogroup_simulation(
        complete_parameters,
        output_abolsute_path,
        threads,
        current_file_path,
        config_file_path,
    )

## to run the get parameters bit
def run_get_parameters(complete_parameters, output_abolsute_path, threads, current_file_path):
    tool = complete_parameters["tool"].lstrip().rstrip().upper()

    # Per-tool validation of optional flags not checked globally.
    needs_proteomes = tool in ("FASTOMA", "BROCCOLI", "SONICPARANOID2", "ORTHOXML")
    needs_tree      = tool in ("BROCCOLI", "SONICPARANOID2", "ORTHOXML")
    missing = []
    if needs_proteomes and not complete_parameters.get("tools_proteomes"):
        missing.append("--tools-proteomes (required for %s)" % tool)
    if needs_tree and not complete_parameters.get("tool_tree"):
        missing.append("--tool-tree (required for %s — tool does not produce a species tree natively)" % tool)
    if missing:
        print("The following flags are missing for tool %s:" % tool)
        for m in missing:
            print("  " + m)
        sys.exit(1)
    get_parameters_dir = os.path.join(
        current_file_path,
        "scripts",
        "GetParameters"
    )
    train_scripts = {
        "FASTOMA": "FastOMATrain.py",
        "OF3": "OrthoFinderTrain.py",
        "BROCCOLI": "BroccoliTrain.py",
        "SONICPARANOID2": "SonicParanoid2Train.py",
        "ORTHOXML": "OrthoXMLTrain.py",
    }
    if tool not in train_scripts:
        print("ERROR: unsupported tool for GetParameters: %s" % tool)
        sys.exit(1)
    train_script = os.path.join(
        get_parameters_dir,
        train_scripts[tool]
    )
    if os.path.isfile(train_script) == False:
        print("ERROR: GetParameters wrapper not found for %s:" % tool)
        print(train_script)
        sys.exit(1)
    print("Running GetParameters for tool: %s" % tool)
    print("Train script: %s" % train_script)
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        train_scripts[tool].replace(".py", ""), train_script
    )
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    mod.run(
        complete_parameters["tool_output"],
        output_abolsute_path,
        threads,
        tools_proteomes=complete_parameters.get("tools_proteomes"),
        tool_tree=complete_parameters.get("tool_tree"),
    )

## to run pfam
def run_pfam(complete_parameters, output_abolsute_path, threads, current_file_path):
    pfam_scripts_dir = os.path.join(
        current_file_path,
        "scripts",
        "PFAM"
    )
    pfam_results_dir = os.path.join(
        current_file_path,
        "PFAM"
    )
    pfam_wrapper = os.path.join(
        pfam_scripts_dir,
        "pfam_wrapper.py"
    )
    if os.path.isfile(pfam_wrapper) == False:
        print("ERROR: PFAM wrapper not found:")
        print(pfam_wrapper)
        sys.exit(1)
    print("Running PFAM workflow via: %s" % pfam_wrapper)
    import importlib.util
    spec = importlib.util.spec_from_file_location("pfam_wrapper", pfam_wrapper)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    mod.run(
        complete_parameters["Proteome"],
        threads,
        pfam_results_dir,
        output_abolsute_path,
    )
        
        
def run_GA(complete_parameters, output_abolsute_path, threads, current_file_path):

    #Saccharomyces_cerevisiae.fa
    #default_species_pick = os.path.join(current_file_path,"PFAM","pfam_genomes","Saccharomyces_cerevisiae")
    #default_pfam_pick =  os.path.join(current_file_path,"pfam_genomes","Saccharomyces_cerevisiae.fa")
    default_species_pick = "Saccharomyces_cerevisiae"
    ### now read pass this data and pathing to the GA function calls...
    ## import GA functionality from scripts pathing..
    import scripts.GeneticAlgorithm.GeneticAlgorithm as GA
    print(GA)
    GA.GA_workflow(complete_parameters, output_abolsute_path, threads, current_file_path,default_species_pick)
    
    
    
    """
    # below is the old RunOrthogroup that didnt use multiproc
    ## main function to simulate an orthogroup
    def RunOrthogroup(outname):
    """
    
def Worflows():
    print("Starting Workflow...")
    args = vars(Programme_Call())
    current_file_path ="/".join(os.path.realpath(__file__).split("/")[:-1])
    function_call,output_abolsute_path,config_file_path,threads,complete_parameters = Process_Args(args,current_file_path)
    print("Running with %s thread(s)" % threads)
    print("This is the data that will be passed to the workflow functions")
    print("Method Call: %s" % function_call)
    print("Output path: %s" % output_abolsute_path)
    print("Config Call: %s" % config_file_path)
    print("Wrapper location: %s" % current_file_path)

    ## if user has called 'orthogroup simulation'
    if function_call == "og_sim":
        run_orthogroup_simulation(
            complete_parameters,
            output_abolsute_path,
            threads,
            current_file_path,
            config_file_path
        )

    ## if user has called 'orthogroup simulation'
    elif function_call == "get_data":
        run_get_parameters(
            complete_parameters,
            output_abolsute_path,
            threads,
            current_file_path
        )

    ## if user has called 'do pfam profiling'
    elif function_call == "PFAM":
            run_pfam(
                complete_parameters,
                output_abolsute_path,
                threads,
                current_file_path
            )
    
    elif function_call == "GA":
        run_GA(complete_parameters,
                        output_abolsute_path,
                        threads,
                        current_file_path)
        
    else:
        print("Issue with function call %s is not a valid function" % function_call)
    
if __name__ == "__main__":
    Worflows()


