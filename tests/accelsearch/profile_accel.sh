#!/usr/bin/env bash
# Time and profile accelsearch (acceleration and jerk searches) on a quiet machine.
#
# Usage:  profile_accel.sh FFTFILE OUTDIR NAME=ACCELSEARCH [NAME=ACCELSEARCH ...]
#
#   FFTFILE  a normalised (rednoised) .fft with its .inf beside it, e.g.
#            Ter5_080912_DM238.50_red.fft
#   OUTDIR   where everything goes (created)
#   NAME=..  the accelsearch binaries to compare (e.g. nocona=build-nocona/src/accelsearch).
#            Each runs with the caller's environment, so run this inside the environment
#            those builds use, and set PRESTO if their wisdom is not installed.
#
# For each binary and search it records: hyperfine timings (1 thread), a flat perf
# profile by symbol and by library, a call-graph profile (to see which accelsearch
# functions call FFTW), -perfinfo (if the binary has it), and md5s of the ACCEL files,
# so the builds can be checked for identical results.  Then thread scaling for the
# first binary.
#
# Environment knobs:
#   SEARCHES="accel8 accel16 jerk"   which searches (see the table below)
#   RUNS=3 WARMUP=1                  hyperfine settings
#   THREADS="1 2 4 8 16"             thread counts for the scaling runs ("" to skip)
#   PERF_FREQ=999                    perf sampling frequency (Hz)
#   XARGS=""                         extra accelsearch arguments for every run
#                                    (e.g. "-flo 100 -fhi 110" for a quick smoke test)
set -euo pipefail

[ $# -ge 3 ] || { sed -n '2,25p' "$0"; exit 1; }
FFT=$(realpath "$1"); INF=${FFT%.fft}.inf; OUT=$(realpath -m "$2"); shift 2
[ -f "$INF" ] || { echo "missing $INF"; exit 1; }
SEARCHES=${SEARCHES:-accel8 accel16 jerk}
RUNS=${RUNS:-3}; WARMUP=${WARMUP:-1}; PERF_FREQ=${PERF_FREQ:-999}
THREADS=${THREADS-1 2 4 8 16}
XARGS=${XARGS:-}
export OPENBLAS_NUM_THREADS=1
declare -A ARGS=(
    [accel8]="-zmax 200 -numharm 8"
    [accel16]="-zmax 200 -numharm 16"
    [jerk]="-zmax 50 -wmax 200 -numharm 8"
)
# A build-tree binary must load its own libpresto, but an activated environment puts
# its lib directory on LD_LIBRARY_PATH, which beats the binary's RUNPATH.  So when a
# libpresto.so sits beside the binary, put that directory first for its runs.
# BIN is the binary, and PRE the "env LD_LIBRARY_PATH=..." prefix (or "env") to run it with.
names=(); declare -A BIN PRE
for nb in "$@"; do
    n=${nb%%=*}; b=$(realpath "${nb#*=}"); names+=("$n"); BIN[$n]=$b
    [ -x "$b" ] || { echo "not executable: $b"; exit 1; }
    PRE[$n]=env
    [ -e "$(dirname "$b")/libpresto.so" ] &&
        PRE[$n]="env LD_LIBRARY_PATH=$(dirname "$b")${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
done
command -v hyperfine >/dev/null || { echo "hyperfine not found"; exit 1; }
command -v perf >/dev/null || { echo "perf not found"; exit 1; }
mkdir -p "$OUT"
base=$(basename "$FFT")

LOG=$OUT/machine.txt
{
    echo "date:      $(date -Is)"
    echo "host:      $(hostname)"
    echo "cpu:       $(lscpu | sed -n 's/^Model name: *//p')"
    echo "governor:  $(cat /sys/devices/system/cpu/cpu0/cpufreq/scaling_governor 2>/dev/null || echo n/a)"
    echo "fft:       $FFT"
    echo "searches:  $SEARCHES"
    echo "runs:      $RUNS (warmup $WARMUP), perf at $PERF_FREQ Hz"
    for n in "${names[@]}"; do
        echo "binary:    $n = ${BIN[$n]}  (run with: ${PRE[$n]})"
        { ${PRE[$n]} ldd "${BIN[$n]}" 2>&1 || true; } | grep -E 'presto|fftw|not a dynamic' |
            sed "s#^#             #" || true
    done
    echo "uptime before: $(uptime)"
} > "$LOG"

rundir() {  # $1 = name, $2 = search: a directory with the data linked in
    local d=$OUT/$1/$2
    mkdir -p "$d"; ln -sf "$FFT" "$INF" "$d/"
    echo "$d"
}

for s in $SEARCHES; do
    a="${ARGS[$s]} $XARGS"
    hf=()
    for n in "${names[@]}"; do
        d=$(rundir "$n" "$s")
        # -perfinfo, if this accelsearch has it (6.3.1 on)
        # (With no arguments it prints its usage and exits non-zero.)
        if { ${PRE[$n]} "${BIN[$n]}" 2>&1 || true; } | grep -q -- -perfinfo; then
            (cd "$d" && ${PRE[$n]} "${BIN[$n]}" -perfinfo $a -flo 100 -fhi 101 "$base" > perfinfo_run.log 2>&1;
                echo "exit status $?" >> perfinfo_run.log) || true
            sed -n '/Performance information/,/FFTW wisdom/p' "$d/perfinfo_run.log" > "$d/perfinfo.txt"
        fi
        hf+=(-n "$n" "cd $d && ${PRE[$n]} ${BIN[$n]} -ncpus 1 $a $base > run.log 2>&1")
    done
    echo "== $s: timing"
    hyperfine --shell=bash --warmup "$WARMUP" --runs "$RUNS" \
        --export-json "$OUT/time_$s.json" --export-markdown "$OUT/time_$s.md" "${hf[@]}" |
        tee "$OUT/time_$s.txt"
    for n in "${names[@]}"; do
        d=$OUT/$n/$s
        md5sum "$d"/*_ACCEL_* 2>/dev/null | sed "s#$d/##" > "$d/accel.md5" || true
        # Flat profile, then a DWARF call graph (no frame pointers in these builds).
        echo "== $s: perf, $n"
        (cd "$d" && ${PRE[$n]} perf record -q -F "$PERF_FREQ" -o perf.data \
            "${BIN[$n]}" -ncpus 1 $a "$base" > perf_run.log 2>&1)
        perf report -i "$d/perf.data" --no-children --sort symbol --stdio \
            --percent-limit 0.3 2>/dev/null | grep -v '^#' | grep -v '^$' > "$d/flat_symbol.txt"
        perf report -i "$d/perf.data" --no-children --sort dso --stdio \
            2>/dev/null | grep -v '^#' | grep -v '^$' > "$d/flat_dso.txt"
        (cd "$d" && ${PRE[$n]} perf record -q -F $((PERF_FREQ / 4)) --call-graph dwarf,16384 \
            -o perf_cg.data "${BIN[$n]}" -ncpus 1 $a "$base" > perf_cg_run.log 2>&1)
        # Which callers the time (inclusive) goes to: the accelsearch-side functions.
        perf report -i "$d/perf_cg.data" --children --sort symbol --stdio \
            --percent-limit 1 -g none 2>/dev/null | grep -v '^#' | grep -v '^$' > "$d/inclusive.txt"
        rm -f "$d/perf_cg.data"
    done
done

# Same candidates from every binary?
{
    echo "ACCEL file check (md5 per search, against ${names[0]}):"
    for s in $SEARCHES; do
        for n in "${names[@]:1}"; do
            if cmp -s "$OUT/${names[0]}/$s/accel.md5" "$OUT/$n/$s/accel.md5"; then r=identical
            else r=DIFFERS; fi
            echo "  $s: $n $r"
        done
    done
} | tee -a "$LOG"

# Thread scaling for the first binary (1 run each after a warmup).
if [ -n "$THREADS" ]; then
    n=${names[0]}
    for s in $SEARCHES; do
        a="${ARGS[$s]} $XARGS"; d=$(rundir "${n}_threads" "$s"); hf=()
        for t in $THREADS; do
            hf+=(-n "$t" "cd $d && ${PRE[$n]} ${BIN[$n]} -ncpus $t $a $base > run_$t.log 2>&1")
        done
        echo "== $s: threads ($n)"
        hyperfine --shell=bash --warmup 1 --runs 2 --export-markdown "$OUT/threads_$s.md" \
            "${hf[@]}" | tee "$OUT/threads_$s.txt"
    done
fi
echo "uptime after:  $(uptime)" >> "$LOG"
echo "results in $OUT"
