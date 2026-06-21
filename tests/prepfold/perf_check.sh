#!/bin/sh
# Performance measurement for prepfold_multi:
#   A) read-once  : one multi pass (3 cands) vs three separate single prepfold
#                   runs -- the multi run pays the 173 MB raw read once, the
#                   singles pay it three times.
#   B) OMP scaling: one 8-candidate multi run at -ncpus 1 / 2 / 4 -- the fold
#                   pass and optimize loop parallelize over candidates, so with
#                   ncand >= ncpus the cores are fully used.
cd "$(dirname "$0")" || exit 2
export PATH="/home/rlynch/sandboxes/presto/build/src:$PATH"
COMMON="-noxwin -noscales -nooffsets -nosearch -fine -n 64 -nsub 64"
DATA=Ter5_080912_short2bits.fits
P=8.66430621957513e-3
PD=-5.01154755640048e-11
now() { date +%s.%N; }
el()  { awk "BEGIN{printf \"%.2f\", $2 - $1}"; }

echo "=== A) read-once: 3-cand multi (ncpus=1) vs 3 single prepfold ==="
t0=$(now)
prepfold_multi $COMMON -ncpus 1 -o perfmulti3 -candfile multi_cands.txt $DATA > perfA_multi.log 2>&1
t1=$(now); multi3=$(el "$t0" "$t1")
s0=$(now)
for dm in 238.30 240.30 242.30; do
    prepfold $COMMON -p $P -pd $PD -dm $dm -o perfA_s$dm $DATA > perfA_s$dm.log 2>&1
done
s1=$(now); singles3=$(el "$s0" "$s1")
echo "  multi(1 pass, ncpus=1) = ${multi3}s"
echo "  3x single (3 passes)   = ${singles3}s"
awk "BEGIN{printf \"  read-once speedup       = %.2fx\n\", $singles3/$multi3}"

echo "=== B) OpenMP scaling: 8-cand multi, ncpus 1/2/4 ==="
prev=""
for nc in 1 2 4; do
    t0=$(now)
    prepfold_multi $COMMON -ncpus $nc -o perf8_$nc -candfile multi_cands_perf.txt $DATA > perfB_nc$nc.log 2>&1
    t1=$(now); e=$(el "$t0" "$t1")
    eval "T$nc=$e"
    echo "  ncpus=$nc : ${e}s"
done
awk "BEGIN{printf \"  speedup ncpus1->2 = %.2fx ; ncpus1->4 = %.2fx\n\", $T1/$T2, $T1/$T4}"

echo "=== B') verify 8-cand parallel run stays correct (ncpus1 vs ncpus4 md5) ==="
ok=1
for cand in perfDM235 perfDM237 perfDM239 perfDM241 perfDM243 perfDM245 perfDM247 perfDM249; do
    m1=$(md5sum perf8_1_${cand}.pfd.bestprof 2>/dev/null | awk '{print $1}')
    m4=$(md5sum perf8_4_${cand}.pfd.bestprof 2>/dev/null | awk '{print $1}')
    [ "$m1" = "$m4" ] && [ -n "$m1" ] || { echo "  MISMATCH: $cand ($m1 vs $m4)"; ok=0; }
done
[ "$ok" -eq 1 ] && echo "  OK: all 8 candidates bit-identical ncpus1 vs ncpus4"
