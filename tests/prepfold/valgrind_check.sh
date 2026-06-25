#!/bin/sh
# Memory + race checks for prepfold_multi's new parallel paths.
#
# Folds only a short slice (-start/-end, small -npart) of the raw data so the
# instrumented runs finish quickly while still exercising the full code path:
# setup -> parallel fold pass (per-thread scratch) -> parallel optimize/output
# (per-candidate struct copies, critical-section writer) -> free.
#
#   1. memcheck, ncpus=1  : leaks + invalid memory, serial path
#   2. memcheck, ncpus=2  : leaks + invalid memory, parallel path (memcheck
#                           serializes threads, so this catches OOB/uninit in
#                           the parallel regions deterministically)
#   3. helgrind, ncpus=2  : data races (note: libgomp produces known false
#                           positives; we grep for races in OUR sources)
cd "$(dirname "$0")" || exit 2
export PATH="/home/rlynch/sandboxes/presto/build/src:$PATH"
DATA=Ter5_080912_short2bits.fits
# short slice keeps instrumented runtime reasonable; still >1 block per part
SLICE="-start 0.0 -end 0.15 -npart 8"
COMMON="-noxwin -noscales -nooffsets -nosearch -fine -n 64 -nsub 64"

echo "### 1. memcheck ncpus=1 (serial) ###"
valgrind --tool=memcheck --leak-check=full --show-leak-kinds=definite,indirect \
    --errors-for-leak-kinds=definite --error-exitcode=99 --log-file=vg_serial.log \
    prepfold_multi $COMMON $SLICE -ncpus 1 -o vgser \
    -candfile multi_cands.txt $DATA > vg_serial.run 2>&1
echo "  memcheck ncpus=1 exit=$?"

echo "### 2. memcheck ncpus=2 (parallel) ###"
valgrind --tool=memcheck --leak-check=full --show-leak-kinds=definite,indirect \
    --errors-for-leak-kinds=definite --error-exitcode=99 --log-file=vg_par.log \
    prepfold_multi $COMMON $SLICE -ncpus 2 -o vgpar \
    -candfile multi_cands.txt $DATA > vg_par.run 2>&1
echo "  memcheck ncpus=2 exit=$?"

echo "### 3. helgrind ncpus=2 (races) ###"
valgrind --tool=helgrind --error-exitcode=98 --log-file=hg_par.log \
    prepfold_multi $COMMON -start 0.0 -end 0.08 -npart 4 -ncpus 2 -o hgpar \
    -candfile multi_cands.txt $DATA > hg_par.run 2>&1
echo "  helgrind ncpus=2 exit=$?"

echo "=== memcheck error/leak summary ==="
for f in vg_serial.log vg_par.log; do
    echo "-- $f --"
    grep -E 'ERROR SUMMARY|definitely lost|indirectly lost' "$f"
done

echo "=== races attributed to prepfold sources (helgrind) ==="
# Data race reports name the racing stack frames; show any that touch our code.
if grep -qE 'prepfold_multi\.c|prepfold_pipeline\.c|fold\.c' hg_par.log; then
    grep -nE 'Possible data race|prepfold_multi\.c|prepfold_pipeline\.c|fold\.c' hg_par.log | head -40
    echo "  ^^ review: races referencing our sources"
else
    echo "  NONE: no helgrind race report references prepfold_multi.c / prepfold_pipeline.c / fold.c"
fi
echo "  (helgrind total race reports: $(grep -c 'Possible data race' hg_par.log))"
