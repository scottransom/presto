#!/bin/sh

# Equivalence + read-once test for prepfold_multi.
#
# The defining correctness property of prepfold_multi is "multi == many singles":
# folding N candidates (each with its own DM) in a SINGLE pass over the raw data
# must reproduce, field for field, what N independent prepfold runs produce.  This
# script proves that for both supported inputs -- raw PSRFITS and PRESTO subbands --
# by folding the candidate files with prepfold_multi and comparing each candidate's
# .bestprof against the matching single prepfold run via compare_bestprof.py.
#
# It also measures wall-clock for the one multi run vs the sum of the single runs,
# turning the structural "the raw data is read once" argument into a measured one
# (the multi run pays the disk read once, the N singles pay it N times).
#
# Run from this directory (it reuses the data downloaded by prepfold_tests.sh and
# creates the subbands itself if they are not already present).

FAILED=0
_SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
LOG=multi_output.txt
: > "$LOG"

# Spin parameters shared by every candidate (must match the *_cands.txt files).
P=8.66430621957513e-3
PD=-5.01154755640048e-11
COMMON="-noxwin -noscales -nooffsets -nosearch -fine -n 64"

# elapsed SECONDS.NS between two `date +%s.%N` stamps, via awk (no bc dependency).
elapsed() { awk "BEGIN{printf \"%.2f\", $2 - $1}"; }

compare() {     # $1 = multi .bestprof, $2 = single .bestprof
    if python3 "${_SCRIPT_DIR}/compare_bestprof.py" "$1" "$2" >> "$LOG" 2>&1; then
        echo "  PASS: $1"
    else
        echo "  FAIL: $1  (vs $2)"
        FAILED=1
    fi
}

# --- Test data (mirror prepfold_tests.sh) ------------------------------------
echo "Getting test data (if needed)."
sh "${_SCRIPT_DIR}/get_test_data.sh" || exit 1
if [ ! -f Ter5_080912_DM242.30.sub00 ]; then
    echo "Making PRESTO subbands (DM 242.3)."
    prepsubband -sub -nsub 32 -subdm 242.3 -nobary -o Ter5_080912 \
        Ter5_080912_short2bits.fits >> "$LOG" 2>&1
fi
DATA=Ter5_080912_short2bits.fits
SUBS="Ter5_080912_DM242.30.sub??"

# --- RAWDATA: prepfold_multi (one pass) vs 3 single prepfold runs -------------
echo "RAWDATA equivalence (PSRFITS, 3 DMs in one pass):"
t0=$(date +%s.%N)
prepfold_multi $COMMON -nsub 64 -candfile "${_SCRIPT_DIR}/multi_cands.txt" \
    $DATA >> "$LOG" 2>&1
t1=$(date +%s.%N)
multi_raw=$(elapsed "$t0" "$t1")

single_raw=0
for pair in 238.30:candDM238 240.30:candDM240 242.30:candDM242; do
    dm=${pair%%:*}; nm=${pair##*:}
    s0=$(date +%s.%N)
    prepfold $COMMON -nsub 64 -p $P -pd $PD -dm $dm -o m_single_$nm \
        $DATA >> "$LOG" 2>&1
    s1=$(date +%s.%N)
    single_raw=$(awk "BEGIN{printf \"%.2f\", $single_raw + ($s1 - $s0)}")
    compare "Ter5_080912_short2bits_${nm}.pfd.bestprof" \
            "m_single_${nm}_8.66ms_Cand.pfd.bestprof"
done
echo "  read-once wall-time: multi=${multi_raw}s  sum-of-singles=${single_raw}s"
awk "BEGIN{ if ($multi_raw < $single_raw) print \"  READ-ONCE OK: multi run was faster than N separate reads\"; else print \"  READ-ONCE WARN: multi was not faster -- check raw I/O is shared\" }"

# --- insubs: prepfold_multi (one pass) vs 3 single prepfold runs --------------
echo "insubs equivalence (PRESTO subbands, 3 residual DMs in one pass):"
prepfold_multi $COMMON -candfile "${_SCRIPT_DIR}/multi_cands_sub.txt" \
    $SUBS >> "$LOG" 2>&1
for pair in 241.30:subDM241 242.30:subDM242 243.30:subDM243; do
    dm=${pair%%:*}; nm=${pair##*:}
    prepfold $COMMON -p $P -pd $PD -dm $dm -o m_single_$nm $SUBS >> "$LOG" 2>&1
    compare "Ter5_080912_DM242.30_${nm}.pfd.bestprof" \
            "m_single_${nm}_8.66ms_Cand.pfd.bestprof"
done

# --- Polyco folds: -timing with par files, and -usepoly -----------------------
# A globular-cluster style timing run: Ter5 A, C, N and O folded with tempo2
# polycos from their par files in one pass (multi_cands_timing.txt has a 'par'
# column), checked against single prepfold -timing runs, by detection, and by
# their TOAs (so get_TOAs.py can time prepfold_multi's polyco folds).
if [ ! -f Ter5_080912_short2bits_rfifind.mask ]; then
    echo "Making the rfifind mask."
    rfifind -noscales -nooffsets -time 2.0 -o Ter5_080912_short2bits \
        $DATA >> "$LOG" 2>&1
fi
TCOMMON="-noxwin -noscales -nooffsets -nsub 64 -mask Ter5_080912_short2bits_rfifind.mask"
echo "Polyco equivalence (-timing, Ter5 A/C/N/O par files in one pass):"
prepfold_multi $TCOMMON -timing -o m_timing -candfile "${_SCRIPT_DIR}/multi_cands_timing.txt" \
    $DATA >> "$LOG" 2>&1
# Only two singles (bright isolated A and the fast binary O): every candidate
# goes through the same code, and the singles are the slow part of the test.
for psr in A O; do
    prepfold $TCOMMON -timing Ter5${psr}.par -o m_single_timing $DATA >> "$LOG" 2>&1
    compare "m_timing_PSR_1748-2446${psr}.pfd.bestprof" \
            "m_single_timing_PSR_1748-2446${psr}.pfd.bestprof"
done
# Detections (reduced chi^2; measured A 1278.4, C 4.00, N 5.01, O 8.99)
for pair in A:1100 C:3.5 N:4.3 O:7.5; do
    psr=${pair%%:*}; min=${pair##*:}
    chi=$(awk '/Reduced chi-sqr/ {print $5}' "m_timing_PSR_1748-2446${psr}.pfd.bestprof")
    if awk "BEGIN{exit !(${chi:-0} > $min)}"; then
        echo "  PASS: Ter5${psr} detected (reduced chi^2 = $chi > $min)"
    else
        echo "  FAIL: Ter5${psr} not detected (reduced chi^2 = $chi, need > $min)"
        FAILED=1
    fi
done
if python3 "${_SCRIPT_DIR}/compare_toas.py" m_timing_PSR_1748-2446O.pfd \
    m_single_timing_PSR_1748-2446O.pfd -t m_timing_PSR_1748-2446O.pfd.bestprof \
    >> "$LOG" 2>&1; then
    echo "  PASS: Ter5O TOAs agree (multi vs single)"
else
    echo "  FAIL: Ter5O TOAs disagree (multi vs single)"
    FAILED=1
fi

# A search candidate folded with polycos made from its barycentric parameters
# (-usepoly; automatic for folds > 2 hr) instead of with bary2topo().
echo "Polyco equivalence (-usepoly):"
printf "# name DM p pd\nupoly 238.30 $P $PD\n" > m_usepoly_cands.txt
prepfold_multi $COMMON -nsub 64 -usepoly -o m_usepoly -candfile m_usepoly_cands.txt \
    $DATA >> "$LOG" 2>&1
prepfold $COMMON -nsub 64 -usepoly -p $P -pd $PD -dm 238.30 -o m_single_usepoly \
    $DATA >> "$LOG" 2>&1
compare m_usepoly_upoly.pfd.bestprof m_single_usepoly_8.66ms_Cand.pfd.bestprof

echo "Finished"
[ "${FAILED}" -eq 0 ] && echo "All prepfold_multi equivalence tests PASSED." \
    || { echo "One or more prepfold_multi equivalence tests FAILED."; exit 1; }

# To clean up temp files, do:
# rm m_single_* m_timing_* m_usepoly* Ter5_080912_short2bits_candDM*.pfd* Ter5_080912_DM242.30_subDM*.pfd* multi_output.txt
