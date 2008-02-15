## Automatically adapted for numpy Apr 14, 2006 by convertcode.py

#
#  From the TEMPO Documentation:
#    
#  The file resid2.tmp contains residuals, etc. in inary format.  
#  Each record contains eight real*8 values:
#    --TOA (MJD, referenced to solar system barycenter)
#    --Postfit residual (pulse phase, from 0 to 1)
#    --Postfit residual (seconds)
#    --Orbital phase (where applicable)
#    --Observing frequency (in barycenter frame)
#    --Weight of point in the fit
#    --Timing uncertainty (according to input file)
#    --Prefit residual (seconds)

import mIO
import numpy as Num

class residuals:
    pass

def read_residuals(filename="resid2.tmp"):
    r = residuals()
    rf = mIO.binary_file(filename)
    # Fortran format: 2 longs (Fortran crap) + 9 doubles = 80 bytes
    r.numTOAs = rf.size()/80
    r.bary_TOA = Num.zeros(r.numTOAs, 'd')
    r.postfit_phs = Num.zeros(r.numTOAs, 'd')
    r.postfit_sec = Num.zeros(r.numTOAs, 'd')
    r.orbit_phs = Num.zeros(r.numTOAs, 'd')
    r.bary_freq = Num.zeros(r.numTOAs, 'd')
    r.weight = Num.zeros(r.numTOAs, 'd')
    r.uncertainty = Num.zeros(r.numTOAs, 'd')
    r.prefit_phs = Num.zeros(r.numTOAs, 'd')
    for ii in range(r.numTOAs):
        struct = rf.fort_read("ddddddddd")
        (r.bary_TOA[ii], 
         r.postfit_phs[ii],
         r.postfit_sec[ii], 
         r.orbit_phs[ii], 
         r.bary_freq[ii], 
         r.weight[ii], 
         r.uncertainty[ii], 
         r.prefit_phs[ii]) = (struct[0], struct[1], struct[2], struct[3], \
                              struct[4], struct[5], struct[6], struct[7])
    rf.close()
    if not Num.nonzero(r.orbit_phs): del r.orbit_phs
    if not Num.nonzero(r.bary_freq): del r.bary_freq
    if not Num.nonzero(r.weight): del r.weight
    r.prefit_sec = r.postfit_sec/r.postfit_phs*r.prefit_phs
    r.uncertainty *= 1.e-6 # Convert uncertainties in usec to sec
    return r

def read_residuals_64bit(filename="resid2.tmp"):
    r = residuals()
    rf = open(filename)
    # Fortran format: 2 8-byte ints (Fortran crap) + 9 doubles = 88 bytes
    rf.seek(0, 2)
    r.numTOAs = rf.tell()/88
    rf.seek(0, 0)
    r.bary_TOA = Num.zeros(r.numTOAs, 'd')
    r.postfit_phs = Num.zeros(r.numTOAs, 'd')
    r.postfit_sec = Num.zeros(r.numTOAs, 'd')
    r.orbit_phs = Num.zeros(r.numTOAs, 'd')
    r.bary_freq = Num.zeros(r.numTOAs, 'd')
    r.weight = Num.zeros(r.numTOAs, 'd')
    r.uncertainty = Num.zeros(r.numTOAs, 'd')
    r.prefit_phs = Num.zeros(r.numTOAs, 'd')
    for ii in range(r.numTOAs):
        block = rf.read(88)
        struct = Num.fromstring(block, 'd')
        (r.bary_TOA[ii], 
         r.postfit_phs[ii],
         r.postfit_sec[ii], 
         r.orbit_phs[ii], 
         r.bary_freq[ii], 
         r.weight[ii], 
         r.uncertainty[ii], 
         r.prefit_phs[ii]) = (struct[1], struct[2], struct[3], struct[4], \
                              struct[5], struct[6], struct[7], struct[8])
    rf.close()
    if not Num.nonzero(r.orbit_phs): del r.orbit_phs
    if not Num.nonzero(r.bary_freq): del r.bary_freq
    if not Num.nonzero(r.weight): del r.weight
    r.prefit_sec = r.postfit_sec/r.postfit_phs*r.prefit_phs
    r.uncertainty *= 1.e-6 # Convert uncertainties in usec to sec
    return r
