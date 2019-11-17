cd /dev/shm
mkdir test
cd test
touch process.stdout
wget http://www.cv.nrao.edu/~sransom/GBT_Lband_PSR.fil
# cp /home/sransom/data/PRESTO_tutorial/GBT_Lband_PSR.fil .
readfile GBT_Lband_PSR.fil 
rfifind -time 1.0 -o Lband GBT_Lband_PSR.fil >> process.stdout
rfifind -nocompute -time 1.0 -freqsig 6.0 -mask Lband_rfifind.mask -o Lband GBT_Lband_PSR.fil >> process.stdout
rfifind_stats.py Lband_rfifind.mask
weights_to_ignorechan.py Lband_rfifind.weights
prepdata -nobary -o Lband_topo_DM0.00 -dm 0.0 -mask Lband_rfifind.mask GBT_Lband_PSR.fil
# exploredat Lband_topo_DM0.00.dat
realfft Lband_topo_DM0.00.dat 
# explorefft Lband_topo_DM0.00.fft
accelsearch -numharm 4 -zmax 0 Lband_topo_DM0.00.dat
cat Lband_topo_DM0.00_ACCEL_0
cp $PRESTO/tests/Lband.birds .
prepdata -o tmp GBT_Lband_PSR.fil | grep Average
DDplan.py -d 500.0 -n 96 -b 96 -t 0.000072 -f 1400.0 -s 32 -r 0.5 -o test_DDplan.ps
ls -l test_DDplan.ps
cp $PRESTO/tests/dedisp.py .
python dedisp.py >> process.stdout 2>&1
mkdir subbands
mv *.sub* subbands/
rm -f tmp.* *.tmp Lband*topo*
ls *.dat | xargs -n 1 realfft >> process.stdout
ls *dat | wc
ls *fft | wc
cp Lband_DM0.00.inf Lband.inf
makezaplist.py Lband.birds 
ls -l *zaplist
ls *.fft | xargs -n 1 zapbirds -zap -zapfile Lband.zaplist -baryv -5.697334e-05 >> process.stdout
ls *.fft | xargs -n 1 accelsearch -zmax 0 -flo 15 >> process.stdout
cp $PRESTO/python/ACCEL_sift.py .
python ACCEL_sift.py
quick_prune_cands.py *DM62*0
quickffdots.py Lband_DM62.00.fft 216.373
gotocand.py -local Lband_DM62.00_ACCEL_0:1 -noxwin
prepfold -noxwin -accelcand 1 -accelfile Lband_DM62.00_ACCEL_0.cand -dm 62 subbands/Lband_DM72.00.sub??
mv subbands/*.pfd* .
prepfold -noxwin -nosearch -n 64 -fine -nsub 96 -p 0.004621638 -dm 62.3 GBT_Lband_PSR.fil 
prepfold -noxwin -nodmsearch -n 64 -fine -nsub 96 -p 0.004621 -o periodsearch -dm 62.3 GBT_Lband_PSR.fil 
pfd_for_timing.py *pfd
single_pulse_search.py *dat >> process.stdout
ls -l Lband_singlepulse.ps
# pygaussfit.py GBT_Lband_PSR_4.62ms_Cand.pfd
# cat > test.gaussians
get_TOAs.py -2 -s 2 -n 16 -d 62.3 -g 0.04 GBT_Lband_PSR_4.62ms_Cand.pfd > xxx.tim
cp $PRESTO/tests/1643-1224.par .
tempo -f 1643-1224.par xxx.tim 
# pyplotres.py 
