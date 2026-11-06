# -*- coding: utf-8 -*-
"""
Created on Wed Jun 17 14:02:26 2026

"""

# -*- coding: utf-8 -*-
"""
Created on Thu Apr 23 11:19:47 2026

@author: JH
"""
import sys
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
import numpy
import pygad
import random
import shutil
import subprocess
import time
import csv
import numpy as np
import scripts.GeneticAlgorithm.ScoreSomeParameters as ScoreSomeParameters
#dup_rate,loss_rate,gbc_val
def write_config_file(complete_parameters,solution,completed_parameters_file):
    input_folder = os.path.abspath(complete_parameters["Orthogroup_train_results"])

    input_parameter_file = param_file = os.path.join(complete_parameters["Orthogroup_train_results"],"TrainingResults","SimulationInputs", "simulation_parameters.txt")


    #ultrametric_tree
    #gap_file
    with open(completed_parameters_file,"a") as comp_para:
        for line in open(input_parameter_file,"r").readlines():
            if line.startswith("max_duplication_rate"):
                comp_para.write("max_duplication_rate=%s\n" % str(solution[0]))
            elif line.startswith("max_loss_rate"):
                comp_para.write("max_loss_rate=%s\n" % str(solution[1]))
            elif line.startswith("gbc"):
                comp_para.write("gbc=%s\n" % str(solution[2]))
            else:
                if not line.startswith("sagephy") and not line.startswith("iqtree"):
                    comp_para.write(line)
        comp_para.write("ultrametric_tree=%s\n" % os.path.join(input_folder,"TrainingResults","SimulationInputs","species_tree.nwk"))
        comp_para.write("gap_file=%s\n" % os.path.join(input_folder,"TrainingResults","SimulationInputs","gap_position_profile_counts.tsv"))


def on_gen(ga):
    percent = round((ga.generations_completed/ga.num_generations)*100)
    print(f"\rCompleted: {percent}% ", end="", flush=True)



def GA_workflow(complete_parameters, output_abolsute_path, threads, current_file_path,default_species_pick):
    ###### import scripts here...
    print("Running GA using %s threads" % str(threads))
    os.mkdir(os.path.join(complete_parameters["output"],"Simulation_Temp_Files"))

    def fitness_func(ga_instance, solution, solution_idx):
        new_solution = [max(0.01, float(x)) for x in solution]
    
        outputdir = None
        try:
            output, outputdir = ScoreSomeParameters.main(
                new_solution[0],
                new_solution[1],
                new_solution[2],
                complete_parameters,
                current_file_path,
                default_species_pick
            )
    
            fitness_results = output.split(",")
            fitness_values = []
            for i in fitness_results:
                I = 1 + float(i)
                fitness_values.append(1 / float(I))
    
            return fitness_values
    
        except Exception as e:
            print("Fitness evaluation failed:", e)
            return [0.0, 0.0, 0.0]
    
        finally:
            if outputdir and os.path.exists(outputdir):
                shutil.rmtree(outputdir, ignore_errors=True)
        #### read function...

    """
    #Example:
    {'output': 'GA_TEST_OUTPUT', 
     'Orthogroups': '50', 
     'generations': '5', 
     'sample_size': '100', 
     'num_parents_mating': '4', 
     'sol_per_pop': '16', 
     'num_genes': '4', 
     'parent_selection_type': 'nsga2', 
     'keep_parents': '1', 
     'crossover_type': 'uniform', 
     'crossover_probability': '0.9', 
     'mutation_type': 'random', 
     'mutation_percent_genes': '10', 
     'save_solutions': '1', 
     'parallel_processing': '16', 
     'Orthogroup_train_results': 'PARAMETERS_FOLDER'}
    
    """
    
    ## deafult..
    function_inputs = [0,0,0] 

    fitness_function = fitness_func
    num_generations = int(complete_parameters["generations"])
    num_parents_mating = int(complete_parameters["num_parents_mating"])
    sol_per_pop = int(complete_parameters["sol_per_pop"])
    num_genes = len(function_inputs)
    #parent_selection_type = "sss"
    parent_selection_type = complete_parameters["parent_selection_type"]
    keep_parents = int(complete_parameters["keep_parents"])
    crossover_type = complete_parameters["crossover_type"]
    keep_elitism=1
    crossover_probability = float(complete_parameters["crossover_probability"])
    mutation_type = complete_parameters["mutation_type"]
    mutation_percent_genes = int(complete_parameters["mutation_percent_genes"])
    sample_size = int(complete_parameters["sample_size"])
    start = time.time()
    time_limit = "time_300"

    ga_instance = pygad.GA(num_generations=num_generations,
                       sample_size = sample_size,
                       num_parents_mating=num_parents_mating,
                       fitness_func=fitness_function,
                       sol_per_pop=sol_per_pop,
                       num_genes=num_genes,
                       parent_selection_type=parent_selection_type,
                       keep_parents=keep_parents,
                       crossover_type=crossover_type,
                       crossover_probability = crossover_probability,
                       mutation_type=mutation_type,
                       keep_elitism = keep_elitism,
                       on_generation=on_gen,
                       mutation_percent_genes=mutation_percent_genes,
                       save_solutions=True,
                       gene_space=[
                            {"low": 0.01, "high": 10.0},
                            {"low": 0.01, "high": 10.0},
                            {"low": 0.01, "high": 10.0},
                        ],
                       save_best_solutions=True,
                       #stop_criteria=time_limit, ## unsupported??
                       suppress_warnings=True,
		       parallel_processing=int(threads)
			)

    ga_instance.run()
    end = time.time()





    report_path = os.path.join(
        complete_parameters["output"],
        "solution_per_generation.csv"
    )



    parameter_names = [
        "parameter_1",
        "parameter_2",
        "parameter_3"
    ]

    fitness_names = [
        "fitness_1",
        "fitness_2",
        "fitness_3"
    ]

    with open(report_path, mode="w", newline="") as csv_file:
        writer = csv.writer(csv_file)

        header = (
            ["generation"]
            + parameter_names
            + fitness_names
            + ["fitness_sum"]
        )

        writer.writerow(header)
        peak_solution = [1,1,1,0]

        for generation, best_solution in enumerate(ga_instance.best_solutions, start=1):

            best_fitness = ga_instance.best_solutions_fitness[generation - 1]

            # Make sure fitness is always list-like.
            # This matters because PyGAD may return scalar fitness for single-objective
            # and list/array fitness for multi-objective.
            if isinstance(best_fitness, (list, tuple, np.ndarray)):
                fitness_values = list(best_fitness)
            else:
                fitness_values = [best_fitness]

            fitness_sum = sum(float(x) for x in fitness_values)

            row = (
                [generation]
                + list(best_solution)
                + fitness_values
                + [fitness_sum]
            )
            if peak_solution[-1] < fitness_sum:
                  peak_solution = [generation,list(best_solution),fitness_values,fitness_sum]
            writer.writerow(row)
    peak_gen,peak_solution,peak_solution_fitness,peak_solution_fitness_sum = peak_solution
    peak_solution = [float(i) for i in peak_solution]
    peak_solution_fitness = [float(i) for i in peak_solution_fitness]
    solution, solution_fitness, solution_idx = ga_instance.best_solution()   
    ## need to find GA peak solution here...	

 
    with open(report_path,"a") as report:
        report.write("*--------- Final Iteration Solution ---------*")
        report.write("Parameters of the last solution : {solution}\n".format(solution=solution))
        report.write("Fitness value of the last solution = {solution_fitness}\n".format(solution_fitness=solution_fitness))
        report.write("*--------- Peak solution of %s at generation %s ---------*" % (str(peak_gen),str(peak_solution_fitness_sum)))
        report.write("Parameters of the peak solution : {peak_solution}\n".format(peak_solution=peak_solution))
        report.write("Fitness value of the peak solution = {peak_solution_fitness}\n".format(peak_solution_fitness=peak_solution_fitness))

        report.write("took : %s\n\n" % str(end- start))
        #report.write("*--------- Pygad Report ---------*")
        #report.write(str(ga_instance.summary()))

    print("\n*------- Final Iteration Solution -------*")
    print("Parameters of the best solution : {solution}".format(solution=solution))
    print("Fitness value of the best solution = {solution_fitness}".format(solution_fitness=solution_fitness))
    print("*--------- Peak solution of %s at generation %s ---------*" % (str(peak_gen),str(peak_solution_fitness_sum)))
    print("Parameters of the peak solution : {peak_solution}\n".format(peak_solution=peak_solution))
    print("Fitness value of the peak solution = {peak_solution_fitness}\n".format(peak_solution_fitness=peak_solution_fitness))

    print("took : %s" % str(end- start))
    print("Saved report to: %s" % report_path)
    completed_parameter_file = os.path.join(complete_parameters['output'],"GA_CompletedParameterFile.txt")
    write_config_file(complete_parameters,solution,completed_parameter_file)

    peak_parameter_file = os.path.join(complete_parameters['output'],"GA_PeakParameterFile.txt")
    write_config_file(complete_parameters,peak_solution,peak_parameter_file)
    print("Wrote Simulation final parameter file to: %s" % peak_parameter_file)
    print("Wrote Simulation peak parameter file to: %s" % completed_parameter_file)



    #print(ga_instance.best_solutions_fitness)
    #return solution, solution_fitness, report_path
    #dup_rate,loss_rate,gbc_val










