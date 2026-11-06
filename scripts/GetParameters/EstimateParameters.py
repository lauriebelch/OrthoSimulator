import os
import random
import subprocess
import statistics
import numpy as np
from multiprocessing import Pool
#from ete3 import Tree
from ete4 import Tree
import math
from itertools import combinations
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent

N_ORTHOGROUPS_TO_SAMPLE = 5000

def write_list(output_dir, filename, values, fmt):
    with open(os.path.join(output_dir, filename), "w") as f:
        for value in values:
            f.write(fmt.format(value) + "\n")

def count_fasta_seqs(filepath):
    count = 0
    with open(filepath, "r") as f:
        for line in f:
            if line.startswith(">"):
                count += 1
    return count


def get_alignment_species_and_gene_counts(fasta_file):
    num_genes = 0
    species_set = set()

    with open(fasta_file) as f:
        for line in f:
            if line.startswith(">"):
                header = line[1:].strip()
                num_genes += 1
                species = header.split("_")[0]
                species_set.add(species)

    return len(species_set), num_genes


def trim_and_summarize(data, lower_percentile=10, upper_percentile=90):
    if not data:
        return None, None, None, None

    data_array = np.array(data)
    lower = np.percentile(data_array, lower_percentile)
    upper = np.percentile(data_array, upper_percentile)
    trimmed = data_array[(data_array >= lower) & (data_array <= upper)]

    trimmed_mean = float(np.mean(trimmed))

    # robust SD using MAD
    median = np.median(trimmed)
    mad = np.median(np.abs(trimmed - median))
    trimmed_sd = float(1.4826 * mad)

    log_mean = math.log(trimmed_mean) if trimmed_mean > 0 else None
    log_sd = (
        math.sqrt(math.log(1 + (trimmed_sd / trimmed_mean) ** 2))
        if trimmed_mean > 0 and trimmed_sd > 0
        else None
    )

    return trimmed_mean, trimmed_sd, log_mean, log_sd


def safe_format(value, digits=4):
    if value is None:
        return "N/A"
    return f"{value:.{digits}f}"


def process_orthogroup(og_file, mean_species_rtt, output_dir, orthogroup_seqs_dir, threads_per_job):
    base = os.path.splitext(og_file)[0]
    print(f"Processing {base}...")

    og_dir = os.path.join(output_dir, "res_files", base)
    os.makedirs(og_dir, exist_ok=True)

    input_path = os.path.join(orthogroup_seqs_dir, og_file)
    aln_file = os.path.join(og_dir, f"{base}.aln.fa")
    trimmed_file = os.path.join(og_dir, f"{base}.trimmed.fa")
    iqtree_prefix = os.path.join(og_dir, base)
    tree_file = f"{iqtree_prefix}.treefile"
    #trimal_stats_file = os.path.join(og_dir, f"{base}_trimal_stats.txt")
    iqtree_file = f"{iqtree_prefix}.iqtree"

    inv_prop = None
    gamma_value = None
    gap_score = None
    median_rtt = None
    median_bl = None
    sigma_tree = None
    wiener = None
    treeness = None
    branch_lengths = []
    num_species = None
    num_genes = None

    try:
        mafft_cmd = [
            "mafft",
            "--maxiterate", "1000",
            "--localpair",
            "--thread", str(threads_per_job),
            "--quiet",
            input_path,
        ]

        with open(aln_file, "w") as out_aln:
            subprocess.run(mafft_cmd, stdout=out_aln, check=True)
        
        trimal_cmd_trim = [
            "trimal",
            "-in", aln_file,
            "-out", trimmed_file,
            "-gappyout",
        ]
        subprocess.run(trimal_cmd_trim, check=True)
        if not os.path.exists(trimmed_file):
            raise FileNotFoundError(f"Trimmed file not found: {trimmed_file}")

        if os.path.getsize(trimmed_file) == 0:
            raise ValueError(f"Trimmed file is empty: {trimmed_file}")

        try:
            num_species, num_genes = get_alignment_species_and_gene_counts(input_path)
        except Exception as e:
            print(f"Error getting species/gene counts for {base}: {e}")

        current_dir = os.path.dirname(os.path.realpath(__file__))
        target_dir = os.path.join(current_dir, '..', '..','bin','iqtree3')
        iqtree_cmd = [
        #    "iqtree3",
            target_dir,
            "-s", trimmed_file,
            "-m", "JTT+G4+I",
            "-pre", iqtree_prefix,
            "-nt", str(threads_per_job),
            "-fast",
            "-quiet",
        ]
        subprocess.run(iqtree_cmd, check=True)

        if os.path.exists(iqtree_file):
            with open(iqtree_file, "r") as f:
                for line in f:
                    if "Number of invariant" in line:
                        try:
                            percent_str = (
                                line.split("(")[-1]
                                .split("%")[-2]
                                .replace("=", "")
                                .strip()
                            )
                            inv_prop = float(percent_str) / 100.0
                        except Exception:
                            continue

                    elif "Gamma shape alpha:" in line:
                        try:
                            gamma_value = float(line.strip().split(":")[-1])
                        except Exception:
                            continue

        if os.path.exists(tree_file):
            try:
                tree = Tree(tree_file, parser=1)
                root = tree.root

                rtt_lengths = [
                    tree.get_distance(root,leaf)
                    for leaf in tree.leaves()
                ]
                median_rtt = statistics.median(rtt_lengths) if rtt_lengths else None

                branch_lengths = [
                    n.dist
                    for n in tree.traverse()
                    if n.dist is not None and n.dist > 0
                ]
                median_bl = statistics.median(branch_lengths) if branch_lengths else None

                ###
                # --- make ultrametric reference tree ---
                tree_ultra = tree.copy()
                tree_ultra.to_ultrametric()
                # map nodes between trees (important: same traversal order)
                nodes = list(tree.traverse())
                nodes_ultra = list(tree_ultra.traverse())
                log_diffs = []

                for n, nu in zip(nodes, nodes_ultra):
                    if n.up is None:
                        continue  # skip root
                    #if n.dist > 0 and n.up.dist > 0 and nu.dist > 0 and nu.up is not None:
                    if n.dist > 0 and nu.dist > 0 and nu.up is not None:
                        if not nu.up.dist:
                            nu.up.dist = 0
                        # observed branch lengths
                        obs = n.dist
                        obs_parent = n.up.dist
                        # expected (ultrametric) branch lengths
                        exp = nu.dist		

                        exp_parent = nu.up.dist if nu.up else None
                        if exp_parent and exp_parent > 0:
                            # relative rates
                            rate = obs / exp
                            rate_parent = obs_parent / exp_parent
                            log_diffs.append(math.log(rate) - math.log(rate_parent))

                sigma_tree = statistics.stdev(log_diffs) if len(log_diffs) > 1 else None
                ###
                leaves = list(tree.leaves())
                wiener = sum(tree.get_distance(a, b) for a, b in combinations(leaves, 2))

                total_bl = sum(branch_lengths)
                internal_bl = sum(
                    n.dist
                    for n in tree.traverse()
                    if not n.is_leaf and n.dist is not None and n.dist > 0
                )

                treeness = internal_bl / total_bl if total_bl > 0 else None

            except Exception as e:
                print(f"Error processing tree file for {base}: {e}")

    except subprocess.CalledProcessError as cpe:
        print(f"Subprocess failed for {base}: {cpe}")

    except Exception as e:
        print(f"Unexpected error processing {base}: {e}")

    return (
        inv_prop,
        gamma_value,
        gap_score,
        median_rtt,
        median_bl,
        sigma_tree,
        wiener,
        treeness,
        branch_lengths,
        num_species,
        num_genes,
    )

### gap stuff

### gap stuff

def read_fasta_alignment(msa_file):
    """
    Read a FASTA alignment into:
        names, sequences
    """
    names = []
    sequences = []

    with open(msa_file) as f:
        name = None
        seq = []

        for line in f:
            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):
                if name is not None:
                    names.append(name)
                    sequences.append("".join(seq))

                name = line[1:]
                seq = []

            else:
                seq.append(line)

        if name is not None:
            names.append(name)
            sequences.append("".join(seq))

    return names, sequences


def get_alignment_gap_profile(msa_file):
    """
    Build one local empirical gap-position profile for one alignment.

    For each alignment column:
        - calculate its relative position along the alignment
        - count how many gap characters occur in that column
        - add that many counts to the corresponding percent bin

    This means a column at 75% along the alignment with 4 gaps
    contributes 4 counts to bin 75 for this alignment only.

    Returns:
        bin_counts:
            raw gap counts per percent bin, indexed 0-99

        probabilities:
            per-bin probabilities normalized within this alignment,
            indexed 0-99
    """
    bin_counts = {
        b: 0
        for b in range(100)
    }

    names, sequences = read_fasta_alignment(msa_file)

    if not sequences:
        probabilities = {
            b: 0.0
            for b in range(100)
        }
        return bin_counts, probabilities

    aln_len = len(sequences[0])

    if aln_len == 0:
        probabilities = {
            b: 0.0
            for b in range(100)
        }
        return bin_counts, probabilities

    for col_idx in range(aln_len):
        gap_count = 0

        for seq in sequences:
            if col_idx < len(seq) and seq[col_idx] == "-":
                gap_count += 1

        if gap_count == 0:
            continue

        percent_position = 100.0 * col_idx / aln_len
        percent_bin = int(math.floor(percent_position))

        if percent_bin > 99:
            percent_bin = 99

        bin_counts[percent_bin] += gap_count

    total_gaps = sum(bin_counts.values())

    if total_gaps > 0:
        probabilities = {
            b: bin_counts[b] / total_gaps
            for b in range(100)
        }
    else:
        probabilities = {
            b: 0.0
            for b in range(100)
        }

    return bin_counts, probabilities


def write_gap_profiles_by_alignment(output_dir, rows):
    """
    Write local empirical gap-position profiles.

    Output:
        gap_position_profiles_by_alignment.tsv

    Columns:
        alignment
        total_gaps
        bin_1
        bin_2
        ...
        bin_100

    Each row is one alignment.

    The bin columns contain probabilities normalized within that alignment.
    Therefore, for each alignment with at least one gap, bin_1 to bin_100
    should sum to approximately 1.0.
    """
    out_file = os.path.join(
        output_dir,
        "gap_position_profile_counts.tsv"
    )

    with open(out_file, "w") as out:
        header = (
            ["alignment", "total_gaps"]
            + [f"bin_{b + 1}" for b in range(100)]
        )

        out.write("\t".join(header) + "\n")

        for alignment_name, total_gaps, probabilities in rows:
            values = [
                alignment_name,
                str(total_gaps),
            ]

            for b in range(100):
                values.append(f"{probabilities.get(b, 0.0):.8f}")

            out.write("\t".join(values) + "\n")

    print(f"Wrote local gap profiles to {out_file}")


def process_folder_gap_profiles(folder, output_dir):
    """
    Process empirical .aln.fa files and create one local
    gap-position profile per alignment.

    Output:
        output_dir/gap_position_profiles_by_alignment.tsv

    Rows:
        one alignment per row

    Columns:
        alignment, total_gaps, bin_1, bin_2, ..., bin_100
    """
    rows = []

    n_files = 0

    for file in Path(folder).rglob("*.aln.fa"):
        bin_counts, probabilities = get_alignment_gap_profile(file)

        alignment_name = file.stem

        # For filenames like OG000001.aln.fa,
        # Path.stem gives OG000001.aln, so remove .aln.
        if alignment_name.endswith(".aln"):
            alignment_name = alignment_name[:-4]

        total_gaps = sum(bin_counts.values())

        rows.append(
            (
                alignment_name,
                total_gaps,
                probabilities,
            )
        )

        n_files += 1

    write_gap_profiles_by_alignment(
        output_dir=output_dir,
        rows=rows,
    )

    print(f"Local gap profiles built from {n_files} alignments")

    return rows

def main(orthofinder_folder, output_folder, n_threads):
    print("Starting parameter estimation...")

    input_base = os.path.abspath(orthofinder_folder)
    output_base = os.path.abspath(output_folder)

    orthogroup_seqs_dir = os.path.join(
        output_base,
        "TrainingResults",
        "Orthogroup_Sequences",
    )

    output_dir = os.path.join(
        output_base,
        "TrainingResults",
        "Orthogroup_analysis",
    )

    duplication_file = os.path.join(
        input_base,
        "Gene_Duplication_Events",
        "Duplications.tsv",
    )

    species_tree_file = os.path.join(
        input_base,
        "Species_Tree",
        "SpeciesTree_rooted_node_labels.txt",
    )

    os.makedirs(output_dir, exist_ok=True)

    threads_per_job = 1
    max_workers = int(n_threads)

    orthogroups = []

    for f in os.listdir(orthogroup_seqs_dir):
        #if f.endswith(".faa"):
        if f.endswith((".fa", ".faa", ".fasta")):
            full_path = os.path.join(orthogroup_seqs_dir, f)
            seq_count = count_fasta_seqs(full_path)

            if 4 <= seq_count < 500:
                orthogroups.append(f)

    sampled = random.sample(
        orthogroups,
        min(N_ORTHOGROUPS_TO_SAMPLE, len(orthogroups)),
    )

    invariable_site_props = []
    gamma_shape_values = []
    gap_scores_weighted = []
    all_branch_lengths = []
    sigma_per_tree = []
    median_rtt_values = []
    wiener_values = []
    treeness_values = []
    num_species_values = []
    num_genes_values = []
    print(species_tree_file)
    species_tree = Tree(species_tree_file, parser=1)
    species_tree.to_ultrametric()
    species_root = species_tree.root
    species_rtt = [
        species_tree.get_distance(leaf,species_root)
        for leaf in species_tree.leaves()
    ]
    mean_species_rtt = statistics.mean(species_rtt)

    with Pool(processes=max_workers) as pool:
        results = pool.starmap(
            process_orthogroup,
            [
                (
                    og,
                    mean_species_rtt,
                    output_dir,
                    orthogroup_seqs_dir,
                    threads_per_job,
                )
                for og in sampled
            ],
        )

    for (
        inv_prop,
        gamma,
        gap_score,
        median_rtt,
        median_bl,
        sigma_tree,
        wiener,
        treeness,
        branch_lengths,
        num_species,
        num_genes,
    ) in results:

        if num_species is not None:
            num_species_values.append(num_species)

        if num_genes is not None:
            num_genes_values.append(num_genes)

        if inv_prop is not None:
            invariable_site_props.append(inv_prop)

        if gamma is not None and gamma <= 10:
            gamma_shape_values.append(gamma)

        if gap_score is not None:
            gap_scores_weighted.append(gap_score)

        if branch_lengths:
            all_branch_lengths.extend(branch_lengths)

        if sigma_tree is not None:
            sigma_per_tree.append(sigma_tree)

        if wiener is not None:
            wiener_values.append(wiener)

        if treeness is not None:
            treeness_values.append(treeness)

        if median_rtt is not None:
            median_rtt_values.append(median_rtt)

    mean_inv, sd_inv, log_mean_inv, log_sd_inv = trim_and_summarize(invariable_site_props)
    mean_gamma, sd_gamma, log_mean_gamma, log_sd_gamma = trim_and_summarize(gamma_shape_values)

    if all_branch_lengths:
        bl_array = np.array(all_branch_lengths)
        lower = np.percentile(bl_array, 5)
        upper = np.percentile(bl_array, 95)
        trimmed_bl = bl_array[(bl_array >= lower) & (bl_array <= upper)]
        mean_bl = np.mean(trimmed_bl)
        start_rate = 2 * mean_bl
    else:
        start_rate = None

    ###
    if sigma_per_tree:
        sigma_array = np.array(sigma_per_tree)
        # remove invalid values
        sigma_array = sigma_array[sigma_array > 0]
        if len(sigma_array) > 0:
            # trim extremes
            lower = np.percentile(sigma_array, 5)
            upper = np.percentile(sigma_array, 95)
            trimmed_sigma = sigma_array[(sigma_array >= lower) & (sigma_array <= upper)]
            # log-transform
            log_sigma = np.log(trimmed_sigma)
            # compute parameters of log-normal
            sigma_log_mean = float(np.mean(log_sigma))
            sigma_log_sd = float(np.std(log_sigma, ddof=1)) if len(log_sigma) > 1 else 0.0
        else:
            sigma_log_mean = None
            sigma_log_sd = None
    else:
        sigma_log_mean = None
        sigma_log_sd = None
    ###

    with open(os.path.join(output_dir, "parameter_estimates.txt"), "w") as f:
        if log_mean_inv is not None:
            f.write(f"Proportion Invariant Sites (Log-transformed) mean: {log_mean_inv:.4f}\n")
            f.write(f"Proportion Invariant Sites (Log-transformed) stddev: {safe_format(log_sd_inv)}\n")
        else:
            f.write("Invariable site proportion: N/A\n")

        if log_mean_gamma is not None:
            f.write(f"Gamma shape (Log-transformed) mean: {log_mean_gamma:.4f}\n")
            f.write(f"Gamma shape (Log-transformed) stddev: {safe_format(log_sd_gamma)}\n")
        else:
            f.write("Gamma shape: N/A\n")

        if start_rate is not None:
            f.write(f"Branch Relaxation Start rate: {start_rate:.6f}\n")
        else:
            f.write("Branch Relaxation Start rate: N/A\n")
        ###
        if sigma_log_mean is not None and sigma_log_sd is not None:
            f.write(f"Branch Relaxation Sigma (log-mean): {sigma_log_mean:.6f}\n")
            f.write(f"Branch Relaxation Sigma (log-sd): {sigma_log_sd:.6f}\n")
        else:
            f.write("Branch Relaxation Sigma: N/A\n")
        ###
    write_list(output_dir, "invariable_site_props.txt", invariable_site_props, "{:.4f}")
    write_list(output_dir, "gamma_shape_values.txt", gamma_shape_values, "{:.4f}")
    write_list(output_dir, "median_rtt.txt", median_rtt_values, "{:.6f}")
    write_list(output_dir, "wiener_index.txt", wiener_values, "{:.6f}")
    write_list(output_dir, "treeness.txt", treeness_values, "{:.6f}")
    write_list(output_dir, "num_species.txt", num_species_values, "{}")
    write_list(output_dir, "num_genes.txt", num_genes_values, "{}")

    if median_rtt_values:
        target_rtt = statistics.median(median_rtt_values)
        scale_factor = target_rtt / mean_species_rtt

        for node in species_tree.traverse():
            if node.dist is not None:
                node.dist *= scale_factor

        species_tree.to_ultrametric()

        output_tree_file = os.path.join(
            output_dir,
            "rescaled_ultrametric_species_tree.nwk",
        )

        species_tree.write(outfile=output_tree_file,parser=1)

    param_file = os.path.join(output_dir, "simulation_parameters.txt")

    with open(param_file, "w") as f:
        f.write(f"prop_invar_mean={safe_format(log_mean_inv)}\n")
        f.write(f"prop_invar_sd={safe_format(log_sd_inv)}\n")
        f.write(f"gamma_shape_mean={safe_format(log_mean_gamma)}\n")
        f.write(f"gamma_shape_sd={safe_format(log_sd_gamma)}\n")

        f.write("max_indel_insert=0.04\n")
        f.write("max_indel_delete=0.04\n")
        f.write("indel_size=1\n")

        f.write("max_duplication_rate=0.5\n")
        f.write("max_loss_rate=0.74\n")
        f.write("max_transfer_rate=0\n")

        f.write("replacement_prob=0.0\n")
        f.write("leaf_sampling_probability=1\n")
        f.write("relax_model=ACRY07\n")

        f.write(f"max_start_rate=1\n")
        ###
        f.write(f"sigma_log_mean={safe_format(sigma_log_mean, 6)}\n")
        f.write(f"sigma_log_sd={safe_format(sigma_log_sd, 6)}\n")
        ###

        f.write("gbc=6.75\n")
        #f.write("sagephy_path=bin/sagephy-1.0.0.jar\n")
        #f.write("iqtree_path=bin/iqtree3\n")

    print("Parameter estimation complete. Summary written to parameter_estimates.txt")
    done_file = os.path.join(output_dir, "DONE")
    with open(done_file, "w") as f:
        f.write("Parameter estimation complete\n")
    
    ### gaps!
    MSA_FOLDER = output_dir

    # -----------------------------
    # PROCESS ONE MSA
    # -----------------------------
    def process_msa_gaps(msa_file):
        """
        Read one FASTA alignment file.
        For each sequence, calculate:
        - number of gap regions
        - median gap-region length
        """
        sequences = []
        # --- read fasta ---
        with open(msa_file) as f:
            name = None
            seq = []
            for line in f:
                line = line.strip()
                if line.startswith(">"):
                    if name is not None:
                        sequences.append("".join(seq))
                    name = line[1:]
                    seq = []
                else:
                    seq.append(line)
            # save final sequence
            if name is not None:
                sequences.append("".join(seq))
        num_gaps_list = []
        median_gap_list = []
        # --- process each sequence ---
        for seq in sequences:
            gap_lengths = []
            current_gap = 0
            gap_start = None
            seq_len = len(seq)
            for i, char in enumerate(seq):
                if char == "-":
                    # start of a new gap region
                    if current_gap == 0:
                        gap_start = i
                    current_gap += 1
                else:
                    # end of a gap region
                    if current_gap > 0:
                        gap_lengths.append(current_gap)
                        current_gap = 0
            # handle a gap region at the end of the sequence
            if current_gap > 0:
                gap_lengths.append(current_gap)
            # only keep sequences that actually have gaps
            if gap_lengths:
                num_gaps_list.append(len(gap_lengths))
                median_gap_list.append(statistics.median(gap_lengths))
        return num_gaps_list, median_gap_list
    
    # -----------------------------
    # PROCESS FOLDER
    # -----------------------------
    def process_folder_gaps(folder):
        """
        Process every .fa alignment file in a folder.
        """
        all_num_gaps = []
        all_median_gaps = []
        for file in Path(folder).rglob("*.aln.fa"):
            n, m = process_msa_gaps(file)
            all_num_gaps.extend(n)
            all_median_gaps.extend(m)
        return all_num_gaps, all_median_gaps
    
    raw_num_gaps, raw_median_gaps = process_folder_gaps(MSA_FOLDER)

    write_list(output_dir, "num_gaps.txt", raw_num_gaps, "{:.4f}")
    write_list(output_dir, "median_gap_size.txt", raw_median_gaps, "{:.4f}")

        # -----------------------------
    # GAP POSITION PROFILES
    # -----------------------------
    process_folder_gap_profiles(
        folder=MSA_FOLDER,
        output_dir=output_dir,
    )

