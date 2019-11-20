from __future__ import print_function
#>>> print fftfit.__doc__
#This module 'fftfit' is auto-generated with f2py (version:2.13.175-1250).
#Functions:
#  zbrent = zbrent(x1,x2,f1,f2,tol,tmp,pha,nsum)
#  dchisqr = dchisqr(tau,tmp,r,nsum)
#  cprof(y,c,amp,pha,nmax=len(y),nh=(len(c)-1))
#  fccf(amp,pha,shift)
#  ffft(d,npts,isign,ireal)
#  fftfit(prof,s,phi,nmax,shift,eshift,snr,esnr,b,errb,ngood)

import numpy as num
from presto.psr_utils import gaussian_profile, TWOPI
from presto.fftfit import cprof, fftfit

template = gaussian_profile(64, 0.5, 0.1)
c,amp,pha = cprof(template)
#pha.savespace()
pha1 = pha[0]
pha = num.fmod(pha-num.arange(1,len(pha)+1)*pha1,TWOPI)
for phs in [0.1, 0.3, 0.7]:
    prof = gaussian_profile(64, phs, 0.1)+num.random.standard_normal(64)
    shift,eshift,snr,esnr,b,errb,ngood = fftfit(prof,amp,pha)
    print("True phs = %f, measured phs = %f +/- %f" % (phs, shift/len(prof),eshift/len(prof)))
