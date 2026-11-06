import os
import argparse
import statistics
from multiprocessing import Pool
from ete4 import Tree

# --- USER SETTINGS ---
MAX_WORKERS = 8  # number of parallel processes (edit as needed)

### python simulate_orthogroup.py -s ../TrainingResults/Orthogroup_analysis/rescaled_ultrametric_species_tree.nwk -f ../UserInput/assign/Saccharomyces_cerevisiae.fasta -PFAM ../UserInput/fungi.pfam.txt -d ../TrainingResults/domain_models.csv -p ../TrainingResults/Orthogroup_analysis/simulation_parameters.txt -o ../sim1 -n 500

# --- FUNCTION: process a single tree ---
def process_tree(fname, tree_dir):
    if not fname.endswith(".tre"):
        return None

    path = os.path.join(tree_dir, fname)

    try:
        tree = Tree(path, parser=1)

        root = tree.root
        rtt = [tree.get_distance(root,leaf) for leaf in tree.leaves()]

        if rtt:
            return statistics.median(rtt)

    except Exception as e:
        print(f"Tree error: {fname} ({e})")

    return None

def GetSimulationStats_main(sim_folder):
#if __name__ == "__main__":

    # --- ARGPARSE ---
    #parser = argparse.ArgumentParser()
    #parser.add_argument("-f", "--folder", required=True, help="Simulation folder")
    #args = parser.parse_args()

    #sim_folder = args.folder

    stats_file = os.path.join(sim_folder, "Simulated_Orthogroup_statistics.txt")
    tree_dir = os.path.join(sim_folder, "tree_files")
    output_dir = os.path.join(sim_folder, "simulation_summaries")

    os.makedirs(output_dir, exist_ok=True)

    # --- 1. PARSE TABLE ---
    num_duplications = []
    num_species = []
    num_genes = []

    with open(stats_file) as f:
        header = f.readline().strip().split("\t")

        required_cols = ["num_duplications", "num_species", "num_genes"]
        for col in required_cols:
            if col not in header:
                raise ValueError(f"Missing column: {col}")

        dup_idx = header.index("num_duplications")
        species_idx = header.index("num_species")
        genes_idx = header.index("num_genes")

        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < len(header):
                continue

            try:
                num_duplications.append(int(parts[dup_idx]))
                num_species.append(int(parts[species_idx]))
                num_genes.append(int(parts[genes_idx]))
            except ValueError:
                continue

    # --- 2. MULTIPROCESS TREE PROCESSING ---
    tree_files = [f for f in os.listdir(tree_dir) if f.endswith(".tre")]

    with Pool(processes=MAX_WORKERS) as pool:
        results = pool.starmap(
            process_tree,
            [(fname, tree_dir) for fname in tree_files]
        )

    # filter out failed results
    median_rtt_values = [r for r in results if r is not None]

    if not median_rtt_values:
        print("Warning: no valid trees processed.")

    # --- 3. WRITE OUTPUT FILES ---
    def write_list(filename, data):
        with open(os.path.join(output_dir, filename), "w") as f:
            for x in data:
                f.write(f"{x}\n")

    write_list("duplication_counts.txt", num_duplications)
    write_list("num_species.txt", num_species)
    write_list("num_genes.txt", num_genes)
    write_list("median_rtt.txt", median_rtt_values)
