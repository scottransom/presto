#!/usr/bin/env bash
# Check that a benchmark will run what you think it will, before spending hours on it.
#
# Usage:  bench_preflight.sh NAME=BINARY [NAME=BINARY ...]
#
#   NAME=..  the PRESTO binaries to be timed (e.g. old=b-old/src/accelsearch), run the
#            way profile_accel.sh runs them: when a libpresto.so sits beside a binary,
#            its directory goes first on LD_LIBRARY_PATH.  Run this inside the same
#            environment (and with the same PRESTO setting) the benchmark will use.
#
# Fails (exit 1) unless:
#   - each binary loads the libpresto beside it (a build tree) or in its own prefix
#     (an install), not an older one that wins on LD_LIBRARY_PATH;
#   - every binary loads the same FFTW library, and that is the expected one:
#     FFTW_EXPECT if set (a path, or part of one, or a version such as 3.3.10),
#     otherwise the one in $CONDA_PREFIX/lib if that exists (the pixi/conda env's);
#   - the FFTW wisdom covers every FFT length each build's accelsearch plans with
#     FFTW_MEASURE (makewisdom -check, from beside each binary when it has -check);
#   - the machine is quiet: fewer than MAXBUSY (default 0.5) CPUs busy over 2 s.
# Linux only (ldd, /proc).  Set PREFLIGHT_WARN_ONLY=1 to report without failing.
set -uo pipefail

[ $# -ge 1 ] || { sed -n '2,22p' "$0"; exit 1; }
MAXBUSY=${MAXBUSY:-0.5}
fails=0
ok()   { echo "  ok    $*"; }
warn() { echo "  WARN  $*"; }
bad()  { echo "  FAIL  $*"; fails=$((fails + 1)); }

# The library a binary resolves for a soname pattern, under its run prefix
resolved() {  # $1 = run prefix, $2 = binary, $3 = soname regex
    $1 ldd "$2" 2>/dev/null | awk -v pat="$3" '$1 ~ pat && $3 ~ /^\// {print $3; exit}'
}
fftw_version() {  # the version string compiled into an FFTW library
    local v
    v=$(strings "$1" 2>/dev/null | grep -E '^fftw-3[.][0-9]+' | head -1)
    echo "${v:-unknown}"
}

if [ -z "${FFTW_EXPECT:-}" ] && [ -n "${CONDA_PREFIX:-}" ] &&
       [ -e "$CONDA_PREFIX/lib/libfftw3f.so.3" ]; then
    FFTW_EXPECT=$(realpath "$CONDA_PREFIX/lib/libfftw3f.so.3")
fi
echo "Preflight on $(hostname -s), $(date -Is)"
echo "  expected FFTW: ${FFTW_EXPECT:-(not set: any, but the same for all binaries)}"
echo "  PRESTO=${PRESTO:-(unset)}"

declare -A FFTWLIB
for nb in "$@"; do
    n=${nb%%=*}; b=$(realpath "${nb#*=}" 2>/dev/null)
    echo "$n: $b"
    [ -x "$b" ] || { bad "not executable"; continue; }
    d=$(dirname "$b"); pre=env
    [ -e "$d/libpresto.so" ] && pre="env LD_LIBRARY_PATH=$d${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"

    # libpresto: beside the binary (build tree) or in the binary's prefix (install)
    lp=$(resolved "$pre" "$b" '^libpresto')
    if [ -z "$lp" ]; then
        bad "libpresto not found"
    else
        lpd=$(dirname "$(realpath "$lp")")
        if [ -e "$d/libpresto.so" ]; then
            [ "$lpd" = "$(realpath "$d")" ] && ok "libpresto $lp" ||
                bad "loads $lp, not the libpresto.so beside it"
        else
            [ "$(dirname "$lpd")" = "$(realpath "$d/..")" ] && ok "libpresto $lp" ||
                bad "loads $lp, from outside its own prefix $(realpath "$d/..")"
        fi
    fi

    # FFTW: the expected one, and the same for every binary
    fl=$(resolved "$pre" "$b" '^libfftw3f[.]')
    if [ -z "$fl" ]; then
        bad "libfftw3f not found"
    else
        FFTWLIB[$n]=$(realpath "$fl"); v=$(fftw_version "$fl")
        if [ -z "${FFTW_EXPECT:-}" ] || [[ ${FFTWLIB[$n]} == *"$FFTW_EXPECT"* ]] ||
               [[ $v == *"$FFTW_EXPECT"* ]]; then
            ok "FFTW $fl ($v)"
        else
            bad "FFTW $fl ($v), expected $FFTW_EXPECT"
        fi
    fi

    # Wisdom, checked by FFTW itself for this build's FFT lengths
    mw=$d/makewisdom
    # (grep -c, not -q: with pipefail, grep -q's early exit can fail the pipeline)
    if [ -x "$mw" ] && [ "$(strings "$mw" | grep -c -- 'makewisdom -check')" -gt 0 ]; then
        out=$($pre "$mw" -check 2>&1); st=$?
        line=$(echo "$out" | grep -E '^(OK|STALE|MISSING|REJECTED)' | head -1)
        wf=$(echo "$out" | sed -n 's/^Wisdom: *//p')
        case $st in
            0) ok "wisdom $wf: $(echo "$line" | sed 's/^OK: *//')" ;;
            *) bad "wisdom $wf: $(echo "$out" | grep -E '^(Missing|STALE|MISSING|REJECTED)' | tr -s ' ' | paste -sd ' ')" ;;
        esac
    else
        warn "no makewisdom with -check beside it: wisdom not checked for this build"
    fi
done

# Every binary on the same FFTW?
if [ ${#FFTWLIB[@]} -gt 1 ] && [ "$(printf '%s\n' "${FFTWLIB[@]}" | sort -u | wc -l)" -gt 1 ]; then
    bad "the binaries load different FFTW libraries:"
    for n in "${!FFTWLIB[@]}"; do echo "          $n: ${FFTWLIB[$n]}"; done
fi

# A quiet machine: busy CPUs over 2 s (iowait, e.g. a RAID resync, counts as idle)
cpustat() { awk '/^cpu /{print $2+$3+$4+$7+$8+$9, $2+$3+$4+$5+$6+$7+$8+$9}' /proc/stat; }
read -r b1 t1 < <(cpustat); sleep 2; read -r b2 t2 < <(cpustat)
busy=$(awk -v b=$((b2 - b1)) -v t=$((t2 - t1)) -v n="$(getconf _NPROCESSORS_ONLN)" \
    'BEGIN {printf "%.2f", (t > 0 ? n * b / t : 0)}')
echo "machine: $(getconf _NPROCESSORS_ONLN) CPUs, load average $(cut -d' ' -f1-3 /proc/loadavg)"
if awk -v x="$busy" -v m="$MAXBUSY" 'BEGIN {exit !(x < m)}'; then
    ok "$busy CPUs busy"
else
    bad "$busy CPUs busy (MAXBUSY=$MAXBUSY).  Busiest now:"
    top -b -n 2 -d 1 -w 200 2>/dev/null | awk '/^top -/{n++} n==2 && $1 ~ /^[0-9]+$/' |
        sort -k9 -nr | head -5 | awk '{printf "          %5s%%  %s (%s)\n", $9, $12, $2}'
fi

if [ $fails -gt 0 ]; then
    echo "Preflight: $fails problem(s)."
    [ "${PREFLIGHT_WARN_ONLY:-0}" = 1 ] && exit 0
    exit 1
fi
echo "Preflight: OK"
