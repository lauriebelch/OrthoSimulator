
import numpy
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
import argparse
from scipy.stats import entropy

# --- ARGPARSE ---





def divergence(dst1, dst2):
    if type(dst1) == list:
        dst1 = numpy.array(dst1)
    if type(dst2) == list:
        dst2 = numpy.array(dst2)

    cmbdst = numpy.concatenate((dst1, dst2))
    Bins = numpy.histogram_bin_edges(cmbdst, bins='fd')

    dist1counts = []
    dist2counts = []

    for b1, b2 in zip(Bins, Bins[1:]):
        mask1dist1 = numpy.abs(dst1) <= b2
        dst1_ = dst1[mask1dist1]
        mask2dist1 = b1 <= numpy.abs(dst1_)
        dst1_res = dst1_[mask2dist1]

        mask1dist2 = numpy.abs(dst2) <= b2
        dst2_ = dst2[mask1dist2]
        mask2dist2 = b1 <= numpy.abs(dst2_)
        dst2_res = dst2_[mask2dist2]

        if len(dst1_res) == 0:
            dist1counts.append(1)
        else:
            dist1counts.append(len(dst1_res))

        if len(dst2_res) == 0:
            dist2counts.append(1)
        else:
            dist2counts.append(len(dst2_res))

    kl_res = entropy(dist1counts, dist2counts, base=2)

    if kl_res == numpy.inf:
        return 1
    else:
        return kl_res

def compare_distributions(emp_folder,sim_folder):
    sim_data_path = os.path.join(sim_folder, "simulation_summaries")
    
    ### empirical distributions
    e_num_genes = numpy.loadtxt(os.path.join(emp_folder, "num_genes.txt"))
    e_num_species = numpy.loadtxt(os.path.join(emp_folder, "num_species.txt"))
    e_median_rtt = numpy.loadtxt(os.path.join(emp_folder, "median_rtt.txt"))
    e_duplications = numpy.loadtxt(os.path.join(emp_folder, "duplications.txt"))
    
    ### simulation distributions
    s_num_genes = numpy.loadtxt(os.path.join(sim_data_path, "num_genes.txt"))
    s_num_species = numpy.loadtxt(os.path.join(sim_data_path, "num_species.txt"))
    s_median_rtt = numpy.loadtxt(os.path.join(sim_data_path, "median_rtt.txt"))
    s_duplications = numpy.loadtxt(os.path.join(sim_data_path, "duplication_counts.txt"))
    
    

    # --- CALCULATE ---
    N_genes = divergence(e_num_genes, s_num_genes)
    N_species = divergence(e_num_species, s_num_species)
    N_rtt = divergence(e_median_rtt, s_median_rtt)
    N_duplications = divergence(e_duplications, s_duplications)
    
    distance = (N_genes + N_species + N_rtt + N_duplications) / 4
    distance = (N_genes + N_species + N_duplications) / 3
    
    score = 1 / (1 + distance)
    return (",".join([str(N_genes),str(N_species),str(N_duplications)]))

    

