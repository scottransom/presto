Binary Optimization Peak Width Results
--------------------------------------

These results are from a Monte Carlo sample of 287 realistic PSR binaries.
Only 112 of these orbits (those in NS-NS systems that had a significant
eccentricity) were used for the Orb e and Orb w width estimation.

The results were:

Logarithmic Fits:

log(Orb p) width =  0.9792168 * log(Orb p/z) - 10.9658871  (chi-sq: 43.697334)
log(Orb x) width =  0.9572412 * log(1.0/z)   +  0.7110553  (chi-sq: 11.204193)
log(Orb t) width =  0.9420009 * log(1.0/z)   -  1.1676730  (chi-sq: 34.879553)

Linear Fits:

Orb p width =  0.0000160 * Orb p/z +  0.0000841  (chi-sq: 0.000816)
Orb x width =  2.2420709 * 1.0/z   +  0.0019139  (chi-sq: 0.004983)
Orb t width =  0.3755975 * 1.0/z   +  0.0002719  (chi-sq: 0.000436)

Note:  The above widths are fractional widths.  Orb t width is the fraction
       of Orb p.  The approximate widths of the Orb e and Orb w peaks are:

Orb e width = 0.016
Orb w width = 0.8 degrees

The raw data is in the file:               testresults.txt
Python scripts for data message are in:    bindata.py
Plots of the above fits are in:            fits.ps and logfits.ps
The Monte Carlo routine itself is in :     montebinopt.py

Scott M. Ransom
5 Nov 1998
