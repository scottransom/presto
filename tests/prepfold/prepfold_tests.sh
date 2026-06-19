#!/bin/sh

# A simple test routine to check the .bestprof files.
# Set STRICT=1 to use exact diff instead of tolerance-based comparison.
STRICT=${STRICT:-0}
FAILED=0

_SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

check_fold() {
    local f="$1"
    local pfd_file="${f%.bestprof}"
    if [ "${STRICT}" = "1" ]; then
        if diff "$f" "goodfolds/$f" >/dev/null 2>&1; then
            echo "PASS (strict): $f"
        else
            echo "FAIL (strict): $f"
            FAILED=1
        fi
    else
        local pfd_arg=""
        [ -f "$pfd_file" ] && pfd_arg="--pfd $pfd_file"
        if python3 "${_SCRIPT_DIR}/compare_bestprof.py" "$f" "goodfolds/$f" $pfd_arg; then
            echo "PASS: $f"
        else
            echo "FAIL: $f"
            FAILED=1
        fi
    fi
}

# Download the file if needed
echo "Getting test data (if needed)."
curl -C - -o Ter5_080912_short2bits.fits https://www.cv.nrao.edu/~sransom/Ter5_080912_short2bits.fits >> output.txt 2>&1

# Dedisperse the data for Ter5N's DM. Do not barycenter.
echo "Processing."
prepdata -nobary -noscales -nooffsets -dm 238.3 -o Ter5_080912_topo_DM238.30 Ter5_080912_short2bits.fits >> output.txt 2>&1

# Fold the topocentric time series using polycos
prepfold -noxwin -timing Ter5N.par Ter5_080912_topo_DM238.30.dat >> output.txt 2>&1
check_fold Ter5_080912_topo_DM238.30_PSR_1748-2446N.pfd.bestprof

# Barycenter the time series.
prepdata -o Ter5_080912_bary_DM238.30 Ter5_080912_topo_DM238.30.dat >> output.txt 2>&1

# Search the time series
accelsearch -zmax 20 Ter5_080912_bary_DM238.30.dat >> output.txt 2>&1

# Fold the time series for Ter5N
prepfold -noxwin -fine -accelcand 2 -accelfile Ter5_080912_bary_DM238.30_ACCEL_20.cand Ter5_080912_bary_DM238.30.dat >> output.txt 2>&1
check_fold Ter5_080912_bary_DM238.30_ACCEL_Cand_2.pfd.bestprof

# Fold the barycentric time series using best-determined spin values
prepfold -noxwin -fine -nosearch -p 8.66430621957513e-3 -pd -5.01154755640048e-11 -n 64 -o Ter5_080912_bary_DM238.30_nosearch Ter5_080912_bary_DM238.30.dat >> output.txt 2>&1
check_fold Ter5_080912_bary_DM238.30_nosearch_8.66ms_Cand.pfd.bestprof

# Fold the barycentric time series but allowing prepfold to optimize p and pdot
prepfold -noxwin -p 8.6642e-3 -pd -3.e-10 -n 64 Ter5_080912_bary_DM238.30.dat >> output.txt 2>&1
check_fold Ter5_080912_bary_DM238.30_8.66ms_Cand.pfd.bestprof

# Make PRESTO subbands at the DM of Ter5A
prepsubband -sub -nsub 32 -subdm 242.3 -nobary -o Ter5_080912 Ter5_080912_short2bits.fits >> output.txt 2>&1

# Fold the PRESTO subbands using polycos
prepfold -noxwin -timing Ter5A.par Ter5_080912_DM242.30.sub?? >> output.txt 2>&1
check_fold Ter5_080912_DM242.30_PSR_1748-2446A.pfd.bestprof

# Run rfifind on the raw data to check masking
rfifind -noscales -nooffsets -time 2.0 -o Ter5_080912_short2bits Ter5_080912_short2bits.fits >> output.txt 2>&1

# Fold the raw data for Ter5A with polycos
prepfold -noxwin -noscales -nooffsets -timing Ter5A.par -mask Ter5_080912_short2bits_rfifind.mask -nsub 64 Ter5_080912_short2bits.fits >> output.txt 2>&1
check_fold Ter5_080912_short2bits_PSR_1748-2446A.pfd.bestprof

# Fold the raw data for Ter5A with a constant period to test barycentering
prepfold -noxwin -noscales -nooffsets -nosearch -fine -p 8.66430621957513e-3 -pd -5.01154755640048e-11 -dm 238.30 -n 64 -nsub 64 Ter5_080912_short2bits.fits >> output.txt 2>&1
check_fold Ter5_080912_short2bits_8.66ms_Cand.pfd.bestprof

# Test folding barycentric events
prepfold -noxwin -fine -nosearch -events -mjds -par M28A.par -n 64 M28A_bary.events >> output.txt 2>&1
check_fold M28A_bary_PSR_1824-2452A.pfd.bestprof

# Stage 3 unit check: the factored read-once/dedisperse-many reader must be
# bit-identical to the original read_subbands().  This needs the non-installed
# 'reader_equiv_check' helper (built via 'meson compile -C build
# reader_equiv_check').  Skipped with a notice if it has not been built.
reader_check_bin=""
for cand in \
    "${PRESTO}/build/src/reader_equiv_check" \
    "${_SCRIPT_DIR}/../../build/src/reader_equiv_check" \
    "$(command -v reader_equiv_check 2>/dev/null)"; do
    [ -n "$cand" ] && [ -x "$cand" ] && { reader_check_bin="$cand"; break; }
done
if [ -n "$reader_check_bin" ]; then
    reader_libdir="$(dirname "$reader_check_bin")"
    for tag in "no mask" "with mask"; do
        if [ "$tag" = "with mask" ]; then
            margs="Ter5_080912_short2bits_rfifind.mask"
        else
            margs=""
        fi
        if LD_LIBRARY_PATH="$reader_libdir:${LD_LIBRARY_PATH}" \
            "$reader_check_bin" Ter5_080912_short2bits.fits $margs >> output.txt 2>&1; then
            echo "PASS: reader_equiv_check ($tag)"
        else
            echo "FAIL: reader_equiv_check ($tag)"
            FAILED=1
        fi
    done
else
    echo "SKIP: reader_equiv_check (build it with 'meson compile -C build reader_equiv_check')"
fi

echo "Finished"
[ "${FAILED}" -eq 0 ] && echo "All tests PASSED." || { echo "One or more tests FAILED."; exit 1; }

# To clean up temp files, do:
# rm *.dat Ter5*inf *rfifind* *.sub* *pfd* *ACCEL* output.txt
