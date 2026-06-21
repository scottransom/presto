#!/bin/sh
# Determinism + -ncpus check for prepfold_multi.
#
# Runs the RAWDATA and insubs candidate folds at ncpus = 1, 4, 8 into distinct
# output roots, then md5-compares each candidate's .bestprof across thread counts.
# Bit-identical across thread counts  ==>  output is deterministic AND -ncpus is
# respected without changing results.  Also tolerance-compares the parallel
# (ncpus=4) output directly against the single-prepfold .bestprof files left by
# prepfold_multi_tests.sh, proving "parallel multi == N singles" directly.
cd "$(dirname "$0")" || exit 2
export PATH="/home/rlynch/sandboxes/presto/build/src:$PATH"
_SCRIPT_DIR="$(pwd)"
COMMON="-noxwin -noscales -nooffsets -nosearch -fine -n 64"
DATA=Ter5_080912_short2bits.fits
SUBS="Ter5_080912_DM242.30.sub??"
FAIL=0

run_set() {   # $1 = label, $2 = root-prefix, $3 = candfile, rest = extra args + data
    label=$1; root=$2; cf=$3; shift 3
    for nc in 1 4 8; do
        prepfold_multi $COMMON -ncpus $nc -o ${root}${nc} \
            -candfile "$cf" "$@" > ${root}_nc${nc}.log 2>&1
        rc=$?
        [ $rc -eq 0 ] || { echo "  ERROR: $label ncpus=$nc exited $rc"; FAIL=1; }
    done
}

triple() {    # $1 = candname, $2 = root-prefix
    f1=${2}1_${1}.pfd.bestprof; f4=${2}4_${1}.pfd.bestprof; f8=${2}8_${1}.pfd.bestprof
    n=$(md5sum "$f1" "$f4" "$f8" 2>/dev/null | awk '{print $1}' | sort -u | wc -l)
    if [ "$n" -eq 1 ]; then
        echo "  DETERMINISTIC ($1): md5 $(md5sum "$f1" | awk '{print $1}') across ncpus 1/4/8"
    else
        echo "  NONDETERMINISTIC ($1):"; md5sum "$f1" "$f4" "$f8"; FAIL=1
    fi
}

echo "=== RAWDATA: ncpus 1/4/8 ==="
run_set RAWDATA detraw multi_cands.txt -nsub 64 $DATA
for cand in candDM238 candDM240 candDM242; do triple $cand detraw; done

echo "=== insubs: ncpus 1/4/8 ==="
run_set insubs detsub multi_cands_sub.txt $SUBS
for cand in subDM241 subDM242 subDM243; do triple $cand detsub; done

echo "=== parallel (ncpus=4) vs single prepfold (direct equivalence) ==="
cmp_one() {   # $1 = parallel bestprof, $2 = single bestprof
    if python3 "${_SCRIPT_DIR}/compare_bestprof.py" "$1" "$2" >/dev/null 2>&1; then
        echo "  PASS: $1 == $2"
    else
        echo "  FAIL: $1 vs $2"; FAIL=1
    fi
}
cmp_one detraw4_candDM238.pfd.bestprof m_single_candDM238_8.66ms_Cand.pfd.bestprof
cmp_one detraw4_candDM240.pfd.bestprof m_single_candDM240_8.66ms_Cand.pfd.bestprof
cmp_one detraw4_candDM242.pfd.bestprof m_single_candDM242_8.66ms_Cand.pfd.bestprof
cmp_one detsub4_subDM241.pfd.bestprof  m_single_subDM241_8.66ms_Cand.pfd.bestprof
cmp_one detsub4_subDM242.pfd.bestprof  m_single_subDM242_8.66ms_Cand.pfd.bestprof
cmp_one detsub4_subDM243.pfd.bestprof  m_single_subDM243_8.66ms_Cand.pfd.bestprof

echo "---"
[ "$FAIL" -eq 0 ] && echo "DETERMINISM+NCPUS: ALL PASS" || echo "DETERMINISM+NCPUS: FAILURES PRESENT"
exit $FAIL
