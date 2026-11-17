#!/usr/bin/env bash
# Time accelsearch across PRESTO releases, for the release-notes speed plot.
#
# Usage:  bench_releases.sh [STAGE ...]      STAGE: build wisdom check run summary
#                                            (default: all five, in that order)
#
# Run it inside the environment whose FFTW and libraries every build should use,
# e.g. on eiger (detached, so it survives the ssh session):
#
#   cd /data1/simd_bench/releases &&
#   (setsid nohup pixi run --manifest-path /data1/environments/pixiPSR/pixi.toml \
#       bash /data1/git/presto/tests/accelsearch/bench_releases.sh > bench.log 2>&1 < /dev/null &)
#
# Each version is built the way a user would build it: the tag's own source (from
# `git archive`, no worktrees), `meson setup` with no options except --prefix (so the
# tag's own default buildtype: release/-O3 up to v5.1.0, debug in v5.2.0,
# debugoptimized/-O2 from v5.3.0; a tag that gcc 14 rejects is rebuilt with only
# its new errors demoted to warnings, see GCC14_COMPAT), installed into its own prefix, and its own
# makewisdom run the way that version documents:
#   - v5.x needs $PRESTO set to its source tree: makewisdom (not installed) is run
#     in $PRESTO/lib, and accelsearch reads $PRESTO/lib/fftw_wisdom.txt and, to
#     label known pulsars in the ACCEL file, $PRESTO/lib/pulsars.cat.  So each v5
#     run has PRESTO=src/<tag>, and only those runs have PRESTO set.
#   - v6.x reads <prefix>/share/presto/fftw_wisdom.txt with PRESTO unset, which is
#     where its installed makewisdom writes.
# Every build compiles and links against the environment's FFTW ($CONDA_PREFIX/lib,
# via PKG_CONFIG_PATH), and every run puts its own prefix's lib dir and then
# $CONDA_PREFIX/lib first on LD_LIBRARY_PATH (v5 installs have no RPATH, and the
# pixi env's LD_LIBRARY_PATH would otherwise load the env's installed libpresto).
# The check stage verifies, per version, with ldd, that accelsearch loads its own
# libpresto and the environment's FFTW and that the wisdom file is there (and, for
# the new version, that `makewisdom -check` passes and -perfinfo agrees).  Every timed
# run's output is also scanned for wisdom warnings.
#
# The new version (NEWREF, labelled NEWLABEL) is built, given wisdom, checked and timed
# first.  It is timed twice: on the *_red.fft file (automatic -nonorm) and on the same
# data under a name without "_red" (block-median normalization, as in older versions);
# its label gets a "-median" suffix for the second.  Older versions use the _red name.
#
# NATIVE_VERSIONS adds "<tag>-native" cases: the same tag built with only
# -Dc_args=-march=native added (its own source tree, prefix and wisdom), to show what
# a user who built for their CPU got.  This matters for v5.3.0-v6.1.1, whose
# accelsearch calls libm's fmaf() in its inner loops unless the compiler may use
# hardware FMA.  ONLY (a regex on case labels) times just some cases, e.g.
#   NATIVE_VERSIONS="v6.1.0 v6.3.3" ONLY='-native$' bench_releases.sh
# appends just the native cases to an existing results.tsv.
#
# Timing: RUNS rounds; each round times every case at every NCPUS, so slow drifts in
# the machine hit all versions alike.  The input is read once beforehand to be in the
# page cache.  The summary takes the fastest wall time per case and thread count.
#
# Environment knobs:
#   BASE=$PWD                     where everything goes (src/ build/ inst/
#                                 run/ logs/, results.tsv, summary.tsv)
#   REPO=<this script's repo>     the git repository the tags come from (needs NEWREF)
#   VERSIONS="v5.0.0 v5.3.1 v6.0.0 v6.1.0 v6.1.2 v6.3.0 v6.3.2"
#                                 released tags (the releases that changed
#                                 accelsearch's internals or build; v5.3.1, not v5.3.0,
#                                 whose rfifind fails to link: no OpenMP in meson)
#   NEWREF=master NEWLABEL=v6.3.3 the upcoming release (any git ref)
#   NATIVE_VERSIONS=""            tags (or NEWLABEL) to also build with -march=native
#   ONLY=""                       regex: time only the cases whose label matches
#   FFTFILE=/data1/cohsearch_tests/Ter5_080912_DM238.50_red.fft   (.inf beside it)
#   ACCELARGS="-zmax 200 -numharm 8"
#   NCPUS="1 4 8"  RUNS=3
#   MAXBUSY=0.5                   passed to bench_preflight.sh's quiet-machine check
set -uo pipefail

SCRIPTDIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
BASE=$(realpath -m "${BASE:-$PWD}")
REPO=${REPO:-$(git -C "$SCRIPTDIR" rev-parse --show-toplevel 2>/dev/null)}
VERSIONS=${VERSIONS:-v5.0.0 v5.3.1 v6.0.0 v6.1.0 v6.1.2 v6.3.0 v6.3.2}
NEWREF=${NEWREF:-master}
NEWLABEL=${NEWLABEL:-v6.3.3}
NATIVE_VERSIONS=${NATIVE_VERSIONS:-}
ONLY=${ONLY:-}
FFTFILE=$(realpath "${FFTFILE:-/data1/cohsearch_tests/Ter5_080912_DM238.50_red.fft}")
ACCELARGS=${ACCELARGS:--zmax 200 -numharm 8}
NCPUS=${NCPUS:-1 4 8}
RUNS=${RUNS:-3}
export MAXBUSY=${MAXBUSY:-0.5}
STAGES=${*:-build wisdom check run summary}
GCC14_COMPAT="-Wno-error=implicit-function-declaration -Wno-error=implicit-int -Wno-error=int-conversion -Wno-error=incompatible-pointer-types"

die() { echo "ERROR: $*" >&2; exit 1; }
log() { echo "[$(date '+%F %T')] $*"; }

[ -n "$REPO" ] || die "set REPO to the PRESTO git repository"
PREFLIGHT=$SCRIPTDIR/bench_preflight.sh
[ -f "$PREFLIGHT" ] || PREFLIGHT=$REPO/tests/accelsearch/bench_preflight.sh
[ -n "${CONDA_PREFIX:-}" ] || die "run this inside the pixi/conda environment (CONDA_PREFIX unset)"
ENVLIB=$CONDA_PREFIX/lib
ENVFFTW=$(realpath "$ENVLIB/libfftw3f.so.3" 2>/dev/null) || die "no $ENVLIB/libfftw3f.so.3"
[ -f "$ENVLIB/pkgconfig/fftw3f.pc" ] || die "no $ENVLIB/pkgconfig/fftw3f.pc"
case $FFTFILE in *_red.fft) ;; *) die "FFTFILE must be a *_red.fft file" ;; esac
[ -f "${FFTFILE%.fft}.inf" ] || die "missing ${FFTFILE%.fft}.inf"
command -v /usr/bin/time > /dev/null || die "needs /usr/bin/time"
NEWSHA=$(git -C "$REPO" rev-parse --verify -q "$NEWREF^{commit}") || die "no ref $NEWREF in $REPO"
mkdir -p "$BASE"/{src,build,inst,run,logs}
cd "$BASE" || exit 1

LABELS="$NEWLABEL $VERSIONS"        # the new version first
for v in $NATIVE_VERSIONS; do LABELS+=" $v-native"; done
tag_of()    { echo "${1%-native}"; }   # a build label's tag (or NEWLABEL)
cargs_of()  { [[ $1 == *-native ]] && echo "-march=native"; }   # its extra c_args
ref_of()    { [ "$(tag_of "$1")" = "$NEWLABEL" ] && echo "$NEWSHA" || tag_of "$1"; }
is_v5()     { [[ $1 == v5.* ]]; }
prefix_of() { echo "$BASE/inst/$1"; }
libdir_of() {  # the lib dir its libpresto was installed into (lib or lib64)
    local f; f=$(ls "$BASE"/inst/"$1"/lib*/libpresto.so 2>/dev/null | head -1)
    [ -n "$f" ] && dirname "$f"
}
wisdom_of() {
    if is_v5 "$1"; then echo "$BASE/src/$1/lib/fftw_wisdom.txt"
    else echo "$(prefix_of "$1")/share/presto/fftw_wisdom.txt"; fi
}
# Run a command in version $1's runtime environment
venv() {
    local v=$1; shift
    local ld; ld=$(libdir_of "$v"):$ENVLIB${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
    if is_v5 "$v"; then
        env LD_LIBRARY_PATH="$ld" PRESTO="$BASE/src/$v" "$@"
    else
        env -u PRESTO LD_LIBRARY_PATH="$ld" "$@"
    fi
}

{
    echo "bench_releases.sh on $(hostname -s), $(date -Is)"
    echo "  stages: $STAGES"
    echo "  versions: $LABELS  ($NEWLABEL = $NEWREF = ${NEWSHA:0:10})"
    echo "  FFTW: $ENVFFTW   compiler: $(cc --version | head -1)"
    echo "  input: $FFTFILE   args: $ACCELARGS   ncpus: $NCPUS   runs: $RUNS"
} | tee -a "$BASE/logs/info.txt"

stage_build() {
    local v ref src bld pre ca
    for v in $LABELS; do
        ref=$(ref_of "$v"); src=$BASE/src/$v; bld=$BASE/build/$v; pre=$(prefix_of "$v")
        ca=$(cargs_of "$v")
        if [ -x "$pre/bin/accelsearch" ]; then log "build $v: already installed"; continue; fi
        log "build $v ($ref${ca:+, c_args $ca})"
        rm -rf "$src" "$bld" "$pre"; mkdir -p "$src"
        git -C "$REPO" archive "$ref" | tar -x -C "$src" || die "git archive $ref failed"
        (
            export PKG_CONFIG_PATH=$ENVLIB/pkgconfig${PKG_CONFIG_PATH:+:$PKG_CONFIG_PATH}
            meson setup "$bld" "$src" --prefix="$pre" ${ca:+-Dc_args="$ca"} &&
                meson compile -C "$bld" && meson install -C "$bld" && exit 0
            # gcc 14 turned some old C warnings into errors (implicit
            # declarations and the like).  None of the default VERSIONS needs
            # it, but an older tag might: demote just those to warnings
            # (diagnostics only, the generated code is unchanged).
            echo "=== bench_releases.sh: retrying with $GCC14_COMPAT"
            meson setup --wipe "$bld" "$src" --prefix="$pre" -Dc_args="${ca:+$ca }$GCC14_COMPAT" &&
                meson compile -C "$bld" && meson install -C "$bld"
        ) > "$BASE/logs/build-$v.log" 2>&1 ||
            die "build of $v failed: see $BASE/logs/build-$v.log"
        grep -q "^=== bench_releases.sh: retrying" "$BASE/logs/build-$v.log" &&
            echo "  (built with -Dc_args=\"$GCC14_COMPAT\")"
        grep -E "Run-time dependency fftw3f|buildtype|optimization" "$bld/meson-logs/meson-log.txt" |
            head -3 | sed "s/^/  /"
        if [ -n "$ca" ]; then
            grep -q -- "$ca" "$bld/build.ninja" || die "$v was not built with $ca"
        fi
        grep -q "$ENVLIB" "$bld/build.ninja" || die "$v did not link the env's FFTW"
    done
}

stage_wisdom() {
    local v w mw
    for v in $LABELS; do
        w=$(wisdom_of "$v")
        if [ -s "$w" ]; then log "wisdom $v: already there ($w)"; continue; fi
        log "wisdom $v -> $w"
        if is_v5 "$v"; then
            # v5: makewisdom is not installed, writes ./fftw_wisdom.txt (FFTW_PATIENT)
            mkdir -p "$(dirname "$w")"
            (cd "$(dirname "$w")" &&
                 venv "$v" "$BASE/build/$v/src/makewisdom") > "$BASE/logs/wisdom-$v.log" 2>&1
        else
            mw=$(prefix_of "$v")/bin/makewisdom
            venv "$v" "$mw" > "$BASE/logs/wisdom-$v.log" 2>&1
        fi
        [ -s "$w" ] || die "makewisdom for $v wrote no $w: see $BASE/logs/wisdom-$v.log"
    done
}

stage_check() {
    local v pre acc lp fl w out fails=0
    log "check: libraries, wisdom, quiet machine"
    for v in $LABELS; do
        pre=$(prefix_of "$v"); acc=$pre/bin/accelsearch
        [ -x "$acc" ] || { echo "  FAIL $v: no $acc"; fails=$((fails + 1)); continue; }
        lp=$(venv "$v" ldd "$acc" | awk '$1 ~ /^libpresto/ {print $3}')
        fl=$(venv "$v" ldd "$acc" | awk '$1 ~ /^libfftw3f[.]/ {print $3}')
        w=$(wisdom_of "$v")
        echo "$v:  libpresto $lp"
        echo "        FFTW      $fl"
        echo "        wisdom    $w"
        [[ $(realpath "$lp" 2>/dev/null) == "$pre"/* ]] ||
            { echo "  FAIL $v: libpresto is not from $pre"; fails=$((fails + 1)); }
        [ "$(realpath "$fl" 2>/dev/null)" = "$ENVFFTW" ] ||
            { echo "  FAIL $v: FFTW is not $ENVFFTW"; fails=$((fails + 1)); }
        [ -s "$w" ] || { echo "  FAIL $v: no wisdom file"; fails=$((fails + 1)); }
        if is_v5 "$v" && [ ! -s "$BASE/src/$v/lib/pulsars.cat" ]; then
            echo "  FAIL $v: no \$PRESTO/lib/pulsars.cat"; fails=$((fails + 1))
        fi
    done
    # The new version can check itself: -perfinfo (libraries, wisdom) and
    # makewisdom -check (the wisdom covers every FFT length it plans), plus
    # bench_preflight.sh's quiet-machine test.
    out=$(venv "$NEWLABEL" "$(prefix_of "$NEWLABEL")/bin/accelsearch" -perfinfo -ncpus 1 2>&1)
    echo "$out" > "$BASE/logs/perfinfo-$NEWLABEL.txt"
    echo "$NEWLABEL -perfinfo:"; echo "$out" | sed 's/^/    /'
    echo "$out" | grep -q "$(prefix_of "$NEWLABEL")/share/presto/fftw_wisdom.txt" ||
        { echo "  FAIL $NEWLABEL: -perfinfo does not show its own wisdom file"; fails=$((fails + 1)); }
    echo "$out" | grep -q "$ENVFFTW\|$ENVLIB/libfftw3f" ||
        { echo "  FAIL $NEWLABEL: -perfinfo does not show the env's FFTW"; fails=$((fails + 1)); }
    venv "$NEWLABEL" bash "$PREFLIGHT" \
         "$NEWLABEL=$(prefix_of "$NEWLABEL")/bin/accelsearch" ||
        { echo "  FAIL: bench_preflight.sh"; fails=$((fails + 1)); }
    [ $fails -eq 0 ] || die "check: $fails problem(s)"
    log "check: OK"
}

# The timed cases: "label version fftname".  The new version also runs on the
# non-_red name, which turns off its automatic -nonorm.
RED=$(basename "$FFTFILE"); NONRED=${RED%_red.fft}.fft
CASES=("$NEWLABEL $NEWLABEL $RED" "$NEWLABEL-median $NEWLABEL $NONRED")
for v in $VERSIONS; do CASES+=("$v $v $RED"); done
for v in $NATIVE_VERSIONS; do CASES+=("$v-native $v-native $RED"); done
if [ -n "$ONLY" ]; then
    sel=(); for c in "${CASES[@]}"; do [[ ${c%% *} =~ $ONLY ]] && sel+=("$c"); done
    [ ${#sel[@]} -gt 0 ] || die "ONLY='$ONLY' matches no case"
    CASES=("${sel[@]}")
fi

stage_run() {
    local res=$BASE/results.tsv c lab v f n r wd tf out st md5 nc
    [ -s "$res" ] || printf "case\tversion\tncpus\trun\twall_s\tuser_s\tsys_s\tmaxrss_kb\tncands\taccel_md5\n" > "$res"
    # Inputs: per-case directories with links (outputs land beside the input)
    for c in "${CASES[@]}"; do
        read -r lab v f <<< "$c"; wd=$BASE/run/$lab; mkdir -p "$wd"
        ln -sf "$FFTFILE" "$wd/$f"; ln -sf "${FFTFILE%.fft}.inf" "$wd/${f%.fft}.inf"
    done
    log "timing: $(for c in "${CASES[@]}"; do printf '%s ' "${c%% *}"; done)"
    log "reading the input into the page cache"; cat "$FFTFILE" > /dev/null
    for r in $(seq 1 "$RUNS"); do
        for c in "${CASES[@]}"; do
            read -r lab v f <<< "$c"; wd=$BASE/run/$lab
            for n in $NCPUS; do
                out=$wd/out_n${n}_r$r.txt; tf=$wd/time_n${n}_r$r.txt
                rm -f "$wd"/*_ACCEL_*
                (cd "$wd" && venv "$v" OMP_NUM_THREADS="$n" /usr/bin/time -f "%e %U %S %M" -o "$tf" \
                     "$(prefix_of "$v")/bin/accelsearch" -ncpus "$n" $ACCELARGS "$f") > "$out" 2>&1
                st=$?
                [ $st -eq 0 ] || die "$lab -ncpus $n run $r failed ($st): see $out"
                if grep -qiE "Couldn.t open.*wisdom|rejected the wisdom|run 'makewisdom'" "$out"; then
                    die "$lab: wisdom problem, see $out"
                fi
                # The new version must use -nonorm on the _red file and block medians otherwise
                if [ "$(tag_of "$v")" = "$NEWLABEL" ]; then
                    if [ "$lab" != "$NEWLABEL-median" ]; then
                        grep -q "as with -nonorm" "$out" || die "$lab did not use -nonorm: $out"
                    else
                        grep -q "as with -nonorm" "$out" && die "$lab used -nonorm: $out"
                    fi
                fi
                md5=$(cat "$wd"/*_ACCEL_*[0-9] 2>/dev/null | md5sum | cut -c1-12)
                nc=$(awk 'NR > 3 && NF == 0 {exit} NR > 3 && $1 ~ /^[0-9]+$/ {n++} END {print n + 0}' \
                         "$wd"/*_ACCEL_*[0-9] 2>/dev/null)
                read -r wall usr sys rss < "$tf"
                printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$lab" "$v" "$n" "$r" \
                       "$wall" "$usr" "$sys" "$rss" "$nc" "$md5" >> "$res"
                log "run $r  $lab  -ncpus $n:  ${wall} s wall, ${usr} s user, $nc cands"
            done
        done
    done
}

stage_summary() {
    local res=$BASE/results.tsv sum=$BASE/summary.tsv
    [ -s "$res" ] || die "no $res"
    # Fastest wall time per case and ncpus, and the speedup over the first
    # released version at the same ncpus
    awk -F'\t' -v base="${VERSIONS%% *}" '
        NR == 1 {next}
        { k = $1 SUBSEP $3; if (!(k in best) || $5 < best[k]) best[k] = $5
          if (!($1 in seen)) {seen[$1] = 1; order[++nc] = $1}
          if (!($3 in seenn)) {seenn[$3] = 1; nn[++nt] = $3}
          cands[k] = $9 }
        END {
          printf "case\tncpus\tbest_wall_s\tspeedup_vs_%s\tncands\n", base
          for (i = 1; i <= nc; i++) for (j = 1; j <= nt; j++) {
              k = order[i] SUBSEP nn[j]; b = base SUBSEP nn[j]
              if (k in best) printf "%s\t%s\t%.2f\t%.2f\t%s\n", order[i], nn[j], best[k],
                                    (b in best ? best[b] / best[k] : 0), cands[k] }
        }' "$res" > "$sum"
    column -t -s$'\t' "$sum"
    echo "(results: $res, summary: $sum; plot with plot_bench_releases.py)"
}

for s in $STAGES; do
    case $s in
        build | wisdom | check | run | summary) "stage_$s" ;;
        *) die "unknown stage $s (build wisdom check run summary)" ;;
    esac
done
log "ALLDONE"
