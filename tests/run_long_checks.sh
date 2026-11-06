#!/usr/bin/env bash
#
# Longer checks for a server run. Run run_example_checks.sh first; this script
# adds bigger and more varied simulations. Everything except the two
# paper_* stages needs only the files in the repository.
#
#   bash tests/run_long_checks.sh --threads 32
#   bash tests/run_long_checks.sh --threads 32 --paper-data /path/to/figure_bundle \
#        --benchmark-data /path/to/manuscript
#
# Stages (all run by default; pick some with --only a,b):
#   scale           one large simulation (--n-large orthogroups), ILS off
#   extremes        simulations at extreme parameter values (high duplication, high loss,
#                   transfers, no indels, ...)
#   species_all     every one of the twelve PFAM species, --n-species orthogroups each
#   ils_range       theta from 0.1 to 1000: lower theta must change the trees more
#   threads         the same small run with 1, 2 and --threads threads
#   paper_repro     repeat the eight simulations of the paper from its settings and compare the
#                   statistics with the paper's data (needs --paper-data)
#   benchmark_repro benchmark.py on the OrthoFinder results of OrthoBench and QfO reproduces the
#                   paper's scores (needs --benchmark-data)
#
# Each stage writes a log to <out>/logs/<stage>.log and the summary is in <out>/summary.txt.

set -u

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT" || exit 1

THREADS=32
N_LARGE=5000
N_SPECIES=200
N_EXTREME=200
N_ILS=500
N_PAPER=3000
OUT=""
ONLY=""
PYTHON=python3
STAGE_TIMEOUT=43200
PAPER_DATA=""
BENCH_DATA=""

usage() {
    sed -n '2,29p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'
    cat <<'EOF'

Options:
  --threads N            threads for every run (default 32)
  --n-large N            orthogroups in the scale stage (default 5000)
  --n-species N          orthogroups per species in species_all (default 200)
  --n-extreme N          orthogroups per run in extremes (default 200)
  --n-ils N              orthogroups per run in ils_range (default 500)
  --n-paper N            orthogroups per run in paper_repro (default 3000)
  --out DIR              output folder (default long_checks_<date and time>)
  --only a,b,c           run only these stages
  --python PATH          python to use (default python3)
  --stage-timeout SECS   time limit for each simulation run (default 43200)
  --paper-data DIR       the figure_bundle folder (SIM_* and Training_* folders)
  --benchmark-data DIR   folder holding OrthoFinder_branch_benchmark and OrthoFinder_branch_benchmark_QFO
  -h, --help             show this help
EOF
}

while [ $# -gt 0 ]; do
    case "$1" in
        --threads) THREADS="$2"; shift 2 ;;
        --n-large) N_LARGE="$2"; shift 2 ;;
        --n-species) N_SPECIES="$2"; shift 2 ;;
        --n-extreme) N_EXTREME="$2"; shift 2 ;;
        --n-ils) N_ILS="$2"; shift 2 ;;
        --n-paper) N_PAPER="$2"; shift 2 ;;
        --out) OUT="$2"; shift 2 ;;
        --only) ONLY="$2"; shift 2 ;;
        --python) PYTHON="$2"; shift 2 ;;
        --stage-timeout) STAGE_TIMEOUT="$2"; shift 2 ;;
        --paper-data) PAPER_DATA="$2"; shift 2 ;;
        --benchmark-data) BENCH_DATA="$2"; shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *) echo "Unknown option: $1"; usage; exit 1 ;;
    esac
done

SIM=OrthoSimulator.py
[ -f "$SIM" ] || { echo "Cannot find OrthoSimulator.py in $REPO_ROOT"; exit 1; }
EXAMPLE_CONFIG=ExampleData/example_simulation_config.txt
[ -n "$OUT" ] || OUT="long_checks_$(date +%Y%m%d_%H%M%S)"
mkdir -p "$OUT/logs" || exit 1
OUT="$(cd "$OUT" && pwd)"
SUMMARY="$OUT/summary.txt"
: > "$SUMMARY"
FAILED=0

##########################################
######### helpers ########################
##########################################

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
    ( "$@" ) > "$log" 2>&1
    local rc=$?
    local secs=$(( $(date +%s) - start ))
    local status=PASS
    if [ $rc -ne 0 ]; then status=FAIL; FAILED=$((FAILED + 1)); fi
    printf "%-5s %-16s %7ss  %s\n" "$status" "$name" "$secs" "$desc" >> "$SUMMARY"
    echo "      $status (${secs}s)   log: $log"
}

skip_stage() {
    local name="$1" why="$2"
    if [ -n "$ONLY" ] && [[ ",$ONLY," != *",$name,"* ]]; then
        return 0
    fi
    printf "%-5s %-16s %7s   %s\n" "SKIP" "$name" "-" "$why" >> "$SUMMARY"
    echo "[$name] skipped: $why"
}

# simulate OUTDIR N CONFIG [extra options]: one run, returns the run's exit status
simulate() {
    local outdir="$1" n="$2" config="$3"
    shift 3
    timeout "$STAGE_TIMEOUT" "$PYTHON" "$SIM" --Orthogroup-simulation --output "$outdir" \
        --Config "$config" --Orthogroups "$n" --threads "$THREADS" "$@"
}

check() { "$PYTHON" tests/check_outputs.py "$@"; }

# every run in a stage is attempted, and the stage fails if any of them failed
STAGE_FAILURES=0
attempt() {
    local label="$1"
    shift
    echo "---- $label"
    if ! "$@"; then
        echo "  FAILED: $label"
        STAGE_FAILURES=$((STAGE_FAILURES + 1))
    fi
}
stage_result() {
    [ "$STAGE_FAILURES" -eq 0 ] || { echo "$STAGE_FAILURES run(s) failed in this stage"; return 1; }
}

##########################################
######### stages #########################
##########################################

sim_and_check() {   # dir n ils theta [extra options]
    local dir="$1" n="$2" ils="$3" theta="$4"
    shift 4
    simulate "$dir" "$n" "$EXAMPLE_CONFIG" "$@" && check sim --dir "$dir" --n "$n" --ils "$ils" --theta "$theta"
}

stage_scale() {
    sim_and_check "$OUT/scale" "$N_LARGE" off 0
    echo "statistics table rows: $(($(wc -l < "$OUT/scale/Simulated_Orthogroup_statistics.txt") - 1))"
}

stage_extremes() {
    STAGE_FAILURES=0
    local n="$N_EXTREME"
    attempt "high duplication"       sim_and_check "$OUT/ex_highdup"   "$n" off 0 --max-duplication-rate 2.0
    attempt "high loss"              sim_and_check "$OUT/ex_highloss"  "$n" off 0 --max-loss-rate 0.95
    attempt "transfers"              sim_and_check "$OUT/ex_transfer"  "$n" off 0 --max-transfer-rate 0.3 --replacement-prob 0.5
    attempt "leaf sampling 0.5"      sim_and_check "$OUT/ex_leaf"      "$n" off 0 --leaf-sampling-probability 0.5
    attempt "no indels"              sim_and_check "$OUT/ex_noindel"   "$n" off 0 --max-indel-insert 0 --max-indel-delete 0
    attempt "frequent indels"        sim_and_check "$OUT/ex_indel"     "$n" off 0 --max-indel-insert 0.3 --max-indel-delete 0.3
    attempt "birth bias 0.5"         sim_and_check "$OUT/ex_gbc_low"   "$n" off 0 --gbc 0.5
    attempt "birth bias 20"          sim_and_check "$OUT/ex_gbc_high"  "$n" off 0 --gbc 20
    attempt "almost no invariant sites" sim_and_check "$OUT/ex_invar"  "$n" off 0 --prop-invar-mean -4 --prop-invar-sd 0.1
    attempt "low gamma shape"        sim_and_check "$OUT/ex_gamma"     "$n" off 0 --gamma-shape-mean -3 --gamma-shape-sd 0.1
    attempt "strong branch-rate variation" sim_and_check "$OUT/ex_sigma" "$n" off 0 --sigma-log-mean 2 --sigma-log-sd 0.5
    stage_result
}

stage_species_all() {
    STAGE_FAILURES=0
    for species in Homo_sapiens Mus_musculus Gallus_gallus Danio_rerio Drosophila_melanogaster \
                   Arabidopsis_thaliana Triticum_aestivum Chlamydomonas_reinhardtii \
                   Saccharomyces_cerevisiae Escherichia_coli Haloferax_volcanii Mnemiopsis_leidyi; do
        attempt "$species" sim_and_check "$OUT/species/$species" "$N_SPECIES" off 0 --species "$species"
    done
    stage_result
}

stage_ils_range() {
    STAGE_FAILURES=0
    local dirs=() thetas=()
    for theta in 0.1 0.3 1 3 10 100 1000; do
        attempt "theta $theta" sim_and_check "$OUT/ils/theta_$theta" "$N_ILS" on "$theta" --ils-theta "$theta"
        dirs+=("$OUT/ils/theta_$theta")
        thetas+=("$theta")
    done
    attempt "tree change falls as theta rises" check theta-trend --dirs "${dirs[@]}" --thetas "${thetas[@]}"
    stage_result
}

stage_threads() {
    STAGE_FAILURES=0
    for t in 1 2 "$THREADS"; do
        attempt "$t thread(s)" bash -c "
            timeout $STAGE_TIMEOUT $PYTHON $SIM --Orthogroup-simulation --output '$OUT/threads/t$t' \
              --Config $EXAMPLE_CONFIG --Orthogroups 50 --threads $t &&
            $PYTHON tests/check_outputs.py sim --dir '$OUT/threads/t$t' --n 50 --ils off"
    done
    stage_result
}

stage_paper_repro() {
    STAGE_FAILURES=0
    for tool in OF3 FastOMA SP2 Broccoli; do
        for dataset in OB QFO; do
            local sim="$PAPER_DATA/SIM_${tool}_${dataset}_10000"
            local training="$PAPER_DATA/Training_${tool}_${dataset}"
            local name="${tool}_${dataset}"
            if [ ! -f "$sim/config.txt" ] || [ ! -d "$training" ]; then
                echo "---- $name: missing $sim or $training"; STAGE_FAILURES=$((STAGE_FAILURES + 1)); continue
            fi
            attempt "$name" bash -c "
                $PYTHON tests/long_check_helpers.py make-config --paper-config '$sim/config.txt' \
                    --training '$training' --n $N_PAPER --out '$OUT/paper/$name.cfg' &&
                timeout $STAGE_TIMEOUT $PYTHON $SIM --Orthogroup-simulation --output '$OUT/paper/$name' \
                    --Config '$OUT/paper/$name.cfg' --threads $THREADS &&
                $PYTHON tests/check_outputs.py sim --dir '$OUT/paper/$name' --n $N_PAPER --ils off &&
                $PYTHON tests/long_check_helpers.py compare-stats \
                    --a '$sim/Simulated_Orthogroup_statistics.txt' --b '$OUT/paper/$name/Simulated_Orthogroup_statistics.txt'"
        done
    done
    stage_result
}

stage_benchmark_repro() {
    STAGE_FAILURES=0
    local bench=benchmarking/benchmark.py
    # scores of OrthoFinder on the OrthoBench and QfO simulations, from the paper's data
    attempt "OrthoBench" bash -c "
        cd '$BENCH_DATA/OrthoFinder_branch_benchmark' &&
        $PYTHON '$REPO_ROOT/$bench' --output-folder reference --tool-output baseline_OF3 --tool orthofinder \
            --csv '$OUT/bench_ob.csv' > /dev/null &&
        cd '$REPO_ROOT' &&
        $PYTHON tests/long_check_helpers.py benchmark-values --csv '$OUT/bench_ob.csv' --expect F1=0.930935 Recall=0.872816"
    attempt "QfO" bash -c "
        cd '$BENCH_DATA/OrthoFinder_branch_benchmark_QFO' &&
        $PYTHON '$REPO_ROOT/$bench' --output-folder reference --tool-output baseline_OF3 --tool orthofinder \
            --csv '$OUT/bench_qfo.csv' > /dev/null &&
        cd '$REPO_ROOT' &&
        $PYTHON tests/long_check_helpers.py benchmark-values --csv '$OUT/bench_qfo.csv' --expect F1=0.9414 Recall=0.9011 --tolerance 5e-4"
    stage_result
}

##########################################
######### run ############################
##########################################

echo "OrthoSimulator long checks"
echo "  repository: $REPO_ROOT"
echo "  output:     $OUT"
echo "  threads:    $THREADS"
echo

run_stage scale           "$N_LARGE orthogroups in one run"                          stage_scale
run_stage extremes        "extreme parameter values ($N_EXTREME orthogroups each)"    stage_extremes
run_stage species_all     "all twelve PFAM species ($N_SPECIES orthogroups each)"     stage_species_all
run_stage ils_range       "ILS from theta 0.1 to 1000 ($N_ILS orthogroups each)"      stage_ils_range
run_stage threads         "the same run with 1, 2 and $THREADS threads"               stage_threads

if [ -n "$PAPER_DATA" ]; then
    run_stage paper_repro "repeat the paper's eight simulations, compare statistics"  stage_paper_repro
else
    skip_stage paper_repro "give --paper-data <figure_bundle folder> to run"
fi

if [ -n "$BENCH_DATA" ]; then
    run_stage benchmark_repro "benchmark.py reproduces the paper's OrthoFinder scores" stage_benchmark_repro
else
    skip_stage benchmark_repro "give --benchmark-data <manuscript folder> to run"
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
