#!/usr/bin/env bash
# Benchmark coherentsearch built with different -march levels, on a quiet machine.
#
# Usage:  bench_march.sh BUILD_DIR FFTFILE [OUTDIR]
#
#   BUILD_DIR  an existing, compiled meson build dir of PRESTO (any c_args; the
#              script strips -march/-mtune from it).  Run this script inside the
#              same environment that build used (e.g. `pixi run bash`), so the
#              same compiler and headers are found.
#   FFTFILE    a normalised .fft with its .inf beside it, e.g.
#              Ter5_080912_DM238.50_red.fft (42 MB).
#   OUTDIR     default ./cs_bench_<host>_<date>
#
# Every variant recompiles only coherentsearch's own sources (the cs_* kernels
# live in the executable, not libpresto) and links against BUILD_DIR's
# libpresto, so libpresto, FFTW and the wisdom are the same for all of them.
#
# Environment knobs:
#   RUNS=5 WARMUP=1   hyperfine settings
#   NCPUS=1           passed as -ncpus
#   VARIANTS="..."    names from the table below (default: all the CPU supports)
#   EXTRA="name=/path/to/coherentsearch ..."  installed binaries to time too
#                     (e.g. a conda-forge one).  They run with the caller's
#                     environment unchanged, so they use their own libs/wisdom.
#   CALLGRIND=1       also count instructions on a 5 Hz slice (needs valgrind)
#   QUICK=1           time the 5 Hz slice instead of the full search (smoke test)
set -euo pipefail

[ $# -ge 2 ] || { sed -n '2,30p' "$0"; exit 1; }
BUILD=$(realpath "$1")
FFT=$(realpath "$2")
INF=${FFT%.fft}.inf
OUT=$(realpath -m "${3:-./cs_bench_$(hostname -s)_$(date +%Y%m%d_%H%M)}")
RUNS=${RUNS:-5}; WARMUP=${WARMUP:-1}; NCPUS=${NCPUS:-1}
SLICE="-lofreq 120 -hifreq 125"
SEARCH=""; [ "${QUICK:-0}" = 1 ] && SEARCH=$SLICE

[ -f "$INF" ] || { echo "missing $INF"; exit 1; }
ls "$BUILD"/src/libpresto.so* >/dev/null 2>&1 || { echo "no libpresto in $BUILD/src: compile it first"; exit 1; }
command -v hyperfine >/dev/null || { echo "hyperfine not found"; exit 1; }

declare -A FLAGS=(
    [nocona]="-march=nocona -mtune=haswell"            # conda-forge's x86-64 flags
    [avx2nofma]="-march=nocona -mtune=haswell -mavx2"  # wider vectors, no FMA
    [v3]="-march=x86-64-v3"                            # AVX2+FMA, portable (Haswell+, Zen+)
    [v4]="-march=x86-64-v4"                            # AVX-512
    [native]="-march=native"
)
if [ -z "${VARIANTS:-}" ]; then
    VARIANTS="nocona avx2nofma v3"
    grep -q avx512f /proc/cpuinfo && VARIANTS="$VARIANTS v4"
    VARIANTS="$VARIANTS native"
fi

mkdir -p "$OUT"
LOG=$OUT/machine.txt
{
    echo "date:      $(date -Is)"
    echo "host:      $(hostname)"
    echo "cpu:       $(lscpu | sed -n 's/^Model name: *//p')"
    echo "governor:  $(cat /sys/devices/system/cpu/cpu0/cpufreq/scaling_governor 2>/dev/null || echo n/a)"
    echo "build:     $BUILD"
    echo "presto:    $(git -C "$BUILD/.." describe --always --dirty 2>/dev/null || echo n/a)"
    echo "fft:       $FFT"
    echo "variants:  $VARIANTS"
    echo "extra:     ${EXTRA:-}"
    echo "runs:      $RUNS (warmup $WARMUP), -ncpus $NCPUS, quick=${QUICK:-0}"
    echo "uptime before: $(uptime)"
} > "$LOG"

# coherentsearch's compile + link commands, as ninja would run them from BUILD.
CMDS=$OUT/cmds.txt
ninja -C "$BUILD" -t commands src/coherentsearch |
    grep -E 'src/coherentsearch\.p/|-o src/coherentsearch ' |
    grep -vE '^(/[^ ]*/)?(meson|python)' > "$CMDS"
CC=$(head -1 "$CMDS" | awk '{print $1}')
echo "compiler:  $($CC --version | head -1)" >> "$LOG"

LDP="$BUILD/src${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
args=()
for v in $VARIANTS; do
    d=$OUT/v_$v; mkdir -p "$d/obj"
    sed -e 's/ -march=[^ ]*//g' -e 's/ -mtune=[^ ]*//g' \
        -e "s#^\([^ ]*\) #\1 ${FLAGS[$v]} #" \
        -e "s#src/coherentsearch\.p/#$d/obj/#g" \
        -e "s#-o src/coherentsearch #-o $d/coherentsearch #" "$CMDS" > "$d/build.sh"
    echo "building $v: ${FLAGS[$v]}"
    (cd "$BUILD" && bash -e "$d/build.sh") > "$d/build.log" 2>&1 ||
        { echo "build of $v failed, see $d/build.log"; exit 1; }
    ln -sf "$FFT" "$INF" "$d/"
    env LD_LIBRARY_PATH="$LDP" ldd "$d/coherentsearch" | grep -E 'presto|fftw' | sed "s#^#  $v: #" >> "$LOG"
    args+=(-n "$v" "cd $d && LD_LIBRARY_PATH=$LDP ./coherentsearch -noprogress -ncpus $NCPUS $SEARCH $(basename "$FFT") > out.txt")
done
for e in ${EXTRA:-}; do
    n=${e%%=*}; b=$(realpath "${e#*=}"); d=$OUT/x_$n; mkdir -p "$d"
    ln -sf "$FFT" "$INF" "$d/"
    ldd "$b" | grep -E 'presto|fftw' | sed "s#^#  $n: #" >> "$LOG"
    args+=(-n "$n" "cd $d && $b -noprogress -ncpus $NCPUS $SEARCH $(basename "$FFT") > out.txt")
done

# One short run each: FFTW must accept the wisdom, or the c2r plans fall back
# to FFTW_ESTIMATE and the FFTs dominate the timing.
for ((i = 1; i < ${#args[@]}; i += 3)); do
    n=${args[i]}; c=${args[i + 1]}
    bash -c "${c/ -noprogress / -noprogress $SLICE }"
    d=$(echo "$c" | sed 's/^cd \([^ ]*\) .*/\1/')
    if grep -qiE 'warning.*(wisdom|makewisdom)|rejected the wisdom' "$d/out.txt"; then
        echo "$n: FFTW wisdom problem:"; cat "$d/out.txt"; exit 1
    fi
done

# Explicit shell: newer hyperfine refuses '&&' without one; bash suits all versions.
hyperfine --shell=bash --warmup "$WARMUP" --runs "$RUNS" \
    --export-json "$OUT/hyperfine.json" --export-markdown "$OUT/hyperfine.md" "${args[@]}" |
    tee "$OUT/hyperfine.txt"
echo "uptime after:  $(uptime)" >> "$LOG"

# Same candidates from every build?  (The full-search .cohout of each.)
ref=$(ls "$OUT"/v_*/*.cohout 2>/dev/null | head -1)
{
    echo "output check against $ref:"
    for f in "$OUT"/v_*/*.cohout "$OUT"/x_*/*.cohout; do
        [ -f "$f" ] || continue
        cmp -s "$ref" "$f" && r=identical || r=DIFFERS
        echo "  $(basename "$(dirname "$f")"): $r"
    done
} | tee -a "$LOG"

if [ "${CALLGRIND:-0}" = 1 ]; then
    echo "instruction counts, $SLICE (total, cs_fill_row, cs_boxcar_batch):" | tee -a "$LOG"
    # valgrind can't decode AVX-512, so skip v4 (and native on AVX-512 CPUs).
    cgv=""
    for v in $VARIANTS; do
        case $v in v4) continue ;; native) grep -q avx512f /proc/cpuinfo && continue ;; esac
        cgv="$cgv $v"
    done
    for v in $cgv; do
        d=$OUT/v_$v
        (cd "$d" && LD_LIBRARY_PATH=$LDP valgrind --tool=callgrind --callgrind-out-file="$d/cg.out" \
            ./coherentsearch -noprogress $SLICE "$(basename "$FFT")" > /dev/null 2>&1) &
    done
    wait
    for v in $cgv; do
        a=$(callgrind_annotate "$OUT/v_$v/cg.out" 2>/dev/null)
        printf "  %-10s %16s %16s %16s\n" "$v" \
            "$(echo "$a" | awk '/PROGRAM TOTALS/{print $1; exit}')" \
            "$(echo "$a" | awk '/cs_fill_row \[/{print $1; exit}')" \
            "$(echo "$a" | awk '/cs_boxcar_batch \[/{print $1; exit}')"
    done | tee -a "$LOG"
fi
echo "results in $OUT"
