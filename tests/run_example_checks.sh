#!/usr/bin/env bash
#
# Checks that run OrthoSimulator end to end on the bundled example data.
#
# Run it from anywhere. It changes to the repository folder itself.
#
#   bash tests/run_example_checks.sh --threads 32
#
# Each stage writes a log to <out>/logs/<stage>.log. A summary table is printed
# at the end and saved to <out>/summary.txt. The exit code is 1 if any stage failed.
#
# Stages that need only the example data (these run by default):
#   environment       Python packages, Java, the bundled IQ-TREE and SaGePhy
#   unit_tests        the unit tests in tests/
#   ils_off           simulation with no --ils-theta: the ILS step must be skipped
#   ils_on            simulation with --ils-theta 10: the ILS step must run
#   theta_sweep       theta = 3, 10, 100: lower theta must change the trees more
#   benchmark_check   benchmark.py scores the true orthogroups against themselves (must be perfect)
#   pfam_species      one small simulation for each of the twelve PFAM species
#
# Optional stages (need data that is not in the repository):
#   baseline_compare  --baseline-ref <branch or commit>  compare with an older version of the code
#   get_parameters    --of3-results, --proteomes, --tool-tree   estimate parameters from OrthoFinder output
#   simulate_trained  runs after get_parameters, simulates from the estimated parameters
#   genetic_algorithm runs after get_parameters, a very short GA run
#   orthoxml          --orthoxml, --orthoxml-proteomes, --orthoxml-tree   estimate parameters from OrthoXML
#   pfam_scan         --with-pfam-scan   the --PFAM workflow (needs pfam_scan.pl, HMMER, Pfam-A and internet)

set -u

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT" || exit 1

THREADS=32
N_OG=100
N_PFAM=10
N_COMPARE=300
OUT=""
ONLY=""
PYTHON=python3
STAGE_TIMEOUT=3600
BASELINE_REF=""
OF3_RESULTS=""
PROTEOMES=""
TOOL_TREE=""
ORTHOXML=""
ORTHOXML_PROTEOMES=""
ORTHOXML_TREE=""
WITH_PFAM_SCAN=0

usage() {
    sed -n '2,30p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'
    cat <<'EOF'

Options:
  --threads N            threads for every run (default 32)
  --n-og N               orthogroups per simulation (default 100)
  --n-pfam N             orthogroups per PFAM species run (default 10)
  --n-compare N          orthogroups per run in the baseline comparison (default 300)
  --out DIR              output folder (default example_checks_<date and time>)
  --only a,b,c           run only these stages
  --python PATH          python to use (default python3)
  --stage-timeout SECS   time limit per simulation run (default 3600)
  --baseline-ref REF     compare with this git branch or commit
  --of3-results DIR      OrthoFinder results folder          (get_parameters)
  --proteomes DIR        proteomes used for that run         (get_parameters)
  --tool-tree FILE       species tree for that run           (get_parameters)
  --orthoxml FILE        an OrthoXML file                    (orthoxml)
  --orthoxml-proteomes DIR
  --orthoxml-tree FILE
  --with-pfam-scan       also run the --PFAM workflow
  -h, --help             show this help
EOF
}

while [ $# -gt 0 ]; do
    case "$1" in
        --threads) THREADS="$2"; shift 2 ;;
        --n-og) N_OG="$2"; shift 2 ;;
        --n-pfam) N_PFAM="$2"; shift 2 ;;
        --n-compare) N_COMPARE="$2"; shift 2 ;;
        --out) OUT="$2"; shift 2 ;;
        --only) ONLY="$2"; shift 2 ;;
        --python) PYTHON="$2"; shift 2 ;;
        --stage-timeout) STAGE_TIMEOUT="$2"; shift 2 ;;
        --baseline-ref) BASELINE_REF="$2"; shift 2 ;;
        --of3-results) OF3_RESULTS="$2"; shift 2 ;;
        --proteomes) PROTEOMES="$2"; shift 2 ;;
        --tool-tree) TOOL_TREE="$2"; shift 2 ;;
        --orthoxml) ORTHOXML="$2"; shift 2 ;;
        --orthoxml-proteomes) ORTHOXML_PROTEOMES="$2"; shift 2 ;;
        --orthoxml-tree) ORTHOXML_TREE="$2"; shift 2 ;;
        --with-pfam-scan) WITH_PFAM_SCAN=1; shift ;;
        -h|--help) usage; exit 0 ;;
        *) echo "Unknown option: $1"; usage; exit 1 ;;
    esac
done

SIM=OrthoSimulator.py
[ -f "$SIM" ] || { echo "Cannot find OrthoSimulator.py in $REPO_ROOT"; exit 1; }

BENCH=""
for candidate in benchmarking/benchmark.py benchmark.py; do
    if [ -f "$candidate" ]; then BENCH="$candidate"; break; fi
done

EXAMPLE_CONFIG=ExampleData/example_simulation_config.txt
[ -n "$OUT" ] || OUT="example_checks_$(date +%Y%m%d_%H%M%S)"
mkdir -p "$OUT/logs" || exit 1
OUT="$(cd "$OUT" && pwd)"
SUMMARY="$OUT/summary.txt"
: > "$SUMMARY"
FAILED=0

##########################################
######### helpers ########################
##########################################

# run a stage: name, description, then the function (and arguments) to run
run_stage() {
    local name="$1" desc="$2"
    shift 2
    if [ -n "$ONLY" ] && [[ ",$ONLY," != *",$name,"* ]]; then
        return 0
    fi
    local log="$OUT/logs/$name.log"
    local start
    start=$(date +%s)
    echo "[$name] $desc"
    ( set -e; "$@" ) > "$log" 2>&1
    local rc=$?
    local secs=$(( $(date +%s) - start ))
    local status=PASS
    if [ $rc -ne 0 ]; then status=FAIL; FAILED=$((FAILED + 1)); fi
    printf "%-5s %-18s %7ss  %s\n" "$status" "$name" "$secs" "$desc" >> "$SUMMARY"
    echo "      $status (${secs}s)   log: $log"
}

skip_stage() {
    local name="$1" why="$2"
    if [ -n "$ONLY" ] && [[ ",$ONLY," != *",$name,"* ]]; then
        return 0
    fi
    printf "%-5s %-18s %7s   %s\n" "SKIP" "$name" "-" "$why" >> "$SUMMARY"
    echo "[$name] skipped: $why"
}

# simulate: output folder, number of orthogroups, then any extra options
simulate() {
    local outdir="$1" n="$2"
    shift 2
    timeout "$STAGE_TIMEOUT" "$PYTHON" "$SIM" --Orthogroup-simulation \
        --output "$outdir" --Config "$EXAMPLE_CONFIG" \
        --Orthogroups "$n" --threads "$THREADS" "$@"
}

check() { "$PYTHON" tests/check_outputs.py "$@"; }

##########################################
######### stages #########################
##########################################

stage_environment() {
    "$PYTHON" --version
    "$PYTHON" - <<'PY'
import importlib, sys
required = ["ete4", "numpy", "scipy", "Bio", "lxml", "mpmath"]
optional = {"pygad": "needed only for the genetic algorithm (--GA)"}
missing = []
for name in required:
    try:
        importlib.import_module(name)
        print("ok       " + name)
    except Exception as error:
        print("MISSING  " + name + "  (" + str(error) + ")")
        missing.append(name)
for name, why in optional.items():
    try:
        importlib.import_module(name)
        print("ok       " + name)
    except Exception:
        print("missing  " + name + "  (" + why + ")")
if missing:
    sys.exit("required packages are missing: " + ", ".join(missing))
PY
    java -version
    bin/iqtree3 --version | head -1
    java -jar bin/sagephy-1.0.0.jar 2>&1 | head -3 || true
    echo "mafft:  $(command -v mafft || echo 'not found (needed by --Get-Parameters)')"
    echo "trimal: $(command -v trimal || echo 'not found (needed by --Get-Parameters)')"
    if command -v git > /dev/null && git rev-parse --git-dir > /dev/null 2>&1; then
        echo "git commit: $(git rev-parse --short HEAD)  branch: $(git rev-parse --abbrev-ref HEAD)"
        echo "uncommitted changes:"; git status --short | head -20
    fi
    echo "threads requested: $THREADS   cores available: $(nproc 2>/dev/null || echo unknown)"
}

stage_unit_tests() {
    "$PYTHON" -B tests/run_all_tests.py
}

stage_ils_off() {
    simulate "$OUT/sim_ils_off" "$N_OG"
    check sim --dir "$OUT/sim_ils_off" --n "$N_OG" --ils off
}

stage_ils_on() {
    simulate "$OUT/sim_ils_on" "$N_OG" --ils-theta 10
    check sim --dir "$OUT/sim_ils_on" --n "$N_OG" --ils on --theta 10
}

stage_theta_sweep() {
    # theta = 1 is left out on purpose: with very strong ILS, large gene families can be slow
    local th
    for th in 3 10 100; do
        simulate "$OUT/sim_theta_$th" "$N_OG" --ils-theta "$th"
    done
    check theta-trend --dirs "$OUT/sim_theta_3" "$OUT/sim_theta_10" "$OUT/sim_theta_100" --thetas 3 10 100
}

stage_benchmark_check() {
    [ -n "$BENCH" ] || { echo "benchmark.py not found"; return 1; }
    [ -f "$OUT/sim_ils_off/Simulated_Orthogroups.txt" ] || { echo "needs the ils_off stage to have run"; return 1; }
    rm -rf "$OUT/perfect_prediction"
    check truth-to-orthofinder --dir "$OUT/sim_ils_off" --out "$OUT/perfect_prediction"
    "$PYTHON" "$BENCH" --output-folder "$OUT/sim_ils_off" --tool-output "$OUT/perfect_prediction" \
        --tool orthofinder --csv "$OUT/perfect_prediction/result.csv"
    check perfect --csv "$OUT/perfect_prediction/result.csv"
}

stage_pfam_species() {
    local species failures=0
    mkdir -p "$OUT/pfam"
    for species in Homo_sapiens Arabidopsis_thaliana Saccharomyces_cerevisiae Escherichia_coli \
                   Mnemiopsis_leidyi Danio_rerio Gallus_gallus Mus_musculus \
                   Chlamydomonas_reinhardtii Drosophila_melanogaster Triticum_aestivum Haloferax_volcanii; do
        echo "=== $species ==="
        if simulate "$OUT/pfam/$species" "$N_PFAM" --species "$species" \
            && check sim --dir "$OUT/pfam/$species" --n "$N_PFAM" --ils off; then
            echo "  OK   $species"
        else
            echo "  FAIL $species"
            failures=$((failures + 1))
        fi
    done
    echo "species that failed: $failures"
    [ "$failures" -eq 0 ]
}

stage_baseline_compare() {
    git rev-parse --git-dir > /dev/null 2>&1 || { echo "not a git repository"; return 1; }
    local dir="$OUT/baseline_code"
    git worktree add --detach "$dir" "$BASELINE_REF"
    local old_sim=OrthoSimulator.py
    local rc=0
    # the runs are not seeded, so both versions get enough orthogroups for the averages to settle
    simulate "$OUT/sim_new_for_compare" "$N_COMPARE"
    ( cd "$dir" && timeout "$STAGE_TIMEOUT" "$PYTHON" "$old_sim" --Orthogroup-simulation \
        --output "$OUT/sim_baseline" --Config "$EXAMPLE_CONFIG" \
        --Orthogroups "$N_COMPARE" --threads "$THREADS" ) || rc=$?
    git worktree remove --force "$dir"
    [ "$rc" -eq 0 ] || return "$rc"
    echo "Comparing this version with ILS off (A) against $BASELINE_REF (B), $N_COMPARE orthogroups each:"
    check compare --a "$OUT/sim_new_for_compare" --b "$OUT/sim_baseline"
}

stage_get_parameters() {
    local run="$OUT/get_parameters"
    "$PYTHON" "$SIM" --Get-Parameters --output "$run" --tool OF3 --tool-output "$OF3_RESULTS" \
        --tool-tree "$TOOL_TREE" --tools-proteomes "$PROTEOMES" --threads "$THREADS"
    local analysis="$run/TrainingResults/Orthogroup_analysis"
    local inputs="$run/TrainingResults/SimulationInputs"
    for f in "$analysis/num_genes.txt" "$analysis/num_species.txt" "$analysis/gamma_shape_values.txt" \
             "$inputs/species_tree.nwk" "$inputs/simulation_parameters.txt" "$inputs/gap_position_profile_counts.tsv"; do
        [ -s "$f" ] || { echo "missing or empty: $f"; return 1; }
        echo "ok  $f"
    done
    echo "orthogroups used to estimate parameters: $(wc -l < "$analysis/num_genes.txt")"
    echo "(the script keeps orthogroups with 4 to 499 genes and samples at most 5000 of them)"
}

stage_simulate_trained() {
    local inputs="$OUT/get_parameters/TrainingResults/SimulationInputs"
    [ -d "$inputs" ] || { echo "needs the get_parameters stage to have run"; return 1; }
    local config="$OUT/trained_config.txt"
    {
        echo "ultrametric_tree=$inputs/species_tree.nwk"
        echo "species=Mnemiopsis_leidyi"
        echo "Orthogroups=$N_OG"
        echo "gap_file=$inputs/gap_position_profile_counts.tsv"
        cat "$inputs/simulation_parameters.txt"
        echo
    } > "$config"
    timeout "$STAGE_TIMEOUT" "$PYTHON" "$SIM" --Orthogroup-simulation --output "$OUT/sim_trained" \
        --Config "$config" --threads "$THREADS"
    check sim --dir "$OUT/sim_trained" --n "$N_OG" --ils off
}

stage_genetic_algorithm() {
    [ -d "$OUT/get_parameters/TrainingResults" ] || { echo "needs the get_parameters stage to have run"; return 1; }
    local config="$OUT/ga_config.txt"
    cat > "$config" <<EOF
output=$OUT/ga_run
Orthogroup_train_results=$OUT/get_parameters
generations=2
sample_size=50
num_parents_mating=4
sol_per_pop=8
num_genes=4
parent_selection_type=nsga2
keep_parents=1
crossover_type=uniform
crossover_probability=0.9
mutation_type=random
mutation_percent_genes=10
save_solutions=1
parallel_processing=$THREADS
Orthogroups=10
EOF
    timeout "$STAGE_TIMEOUT" "$PYTHON" "$SIM" --GA --Config "$config" --threads "$THREADS"
    echo "files written by the GA:"
    ls -R "$OUT/ga_run" | head -40
}

stage_orthoxml() {
    "$PYTHON" "$SIM" --Get-Parameters --output "$OUT/get_parameters_orthoxml" --tool ORTHOXML \
        --tool-output "$ORTHOXML" --tools-proteomes "$ORTHOXML_PROTEOMES" --tool-tree "$ORTHOXML_TREE" \
        --threads "$THREADS"
    local analysis="$OUT/get_parameters_orthoxml/TrainingResults/Orthogroup_analysis"
    [ -s "$analysis/num_genes.txt" ] || { echo "missing or empty: $analysis/num_genes.txt"; return 1; }
    echo "orthogroups used to estimate parameters: $(wc -l < "$analysis/num_genes.txt")"
    echo "duplications recorded (must not all be zero): $(awk '{s+=$1} END {print s}' "$analysis/duplications.txt")"
}

stage_pfam_scan() {
    timeout "$STAGE_TIMEOUT" "$PYTHON" "$SIM" --PFAM --output "$OUT/pfam_run" \
        --Proteome ExampleData/pfam_genome.fa --threads "$THREADS"
    ls -R "$OUT/pfam_run" | head -40
}

##########################################
######### run ############################
##########################################

echo "OrthoSimulator example checks"
echo "  repository: $REPO_ROOT"
echo "  output:     $OUT"
echo "  threads:    $THREADS    orthogroups per run: $N_OG"
echo

run_stage environment      "Python packages, Java, IQ-TREE, SaGePhy"                stage_environment
run_stage unit_tests       "unit tests"                                            stage_unit_tests
run_stage ils_off          "simulation with ILS off (the default)"                  stage_ils_off
run_stage ils_on           "simulation with ILS on (theta = 10)"                    stage_ils_on
run_stage theta_sweep      "ILS strength follows theta (3, 10, 100)"                stage_theta_sweep
run_stage benchmark_check  "benchmark.py gives a perfect score to the true answer"  stage_benchmark_check
run_stage pfam_species     "all twelve PFAM species simulate"                       stage_pfam_species

if [ -n "$BASELINE_REF" ]; then
    run_stage baseline_compare "compare with $BASELINE_REF"                         stage_baseline_compare
else
    skip_stage baseline_compare "give --baseline-ref <branch or commit> to run"
fi

if [ -n "$OF3_RESULTS" ] && [ -n "$PROTEOMES" ] && [ -n "$TOOL_TREE" ]; then
    run_stage get_parameters   "estimate parameters from OrthoFinder output"        stage_get_parameters
    run_stage simulate_trained "simulate from the estimated parameters"             stage_simulate_trained
    run_stage genetic_algorithm "short genetic algorithm run"                       stage_genetic_algorithm
else
    skip_stage get_parameters    "give --of3-results, --proteomes and --tool-tree to run"
    skip_stage simulate_trained  "needs get_parameters"
    skip_stage genetic_algorithm "needs get_parameters"
fi

if [ -n "$ORTHOXML" ] && [ -n "$ORTHOXML_PROTEOMES" ] && [ -n "$ORTHOXML_TREE" ]; then
    run_stage orthoxml "estimate parameters from an OrthoXML file"                  stage_orthoxml
else
    skip_stage orthoxml "give --orthoxml, --orthoxml-proteomes and --orthoxml-tree to run"
fi

if [ "$WITH_PFAM_SCAN" -eq 1 ]; then
    run_stage pfam_scan "the --PFAM workflow on the example proteome"               stage_pfam_scan
else
    skip_stage pfam_scan "give --with-pfam-scan to run"
fi

echo
echo "================ summary ================"
cat "$SUMMARY"
echo "========================================="
if [ "$FAILED" -eq 0 ]; then
    echo "All stages that ran passed."
else
    echo "$FAILED stage(s) FAILED. Read the logs in $OUT/logs/"
fi
echo "Summary: $SUMMARY. Logs are in the logs folder of the output."
[ "$FAILED" -eq 0 ]
