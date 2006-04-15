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
import numpy as Numeric

class residuals:
    pass

def read_residuals():
    r = residuals()
    rf = mIO.binary_file("resid2.tmp")
    # Fortran format: 2 longs (Fortran crap) + 9 doubles = 80 bytes
    r.numTOAs = rf.size()/88
    r.bary_TOA = Numeric.zeros(r.numTOAs, 'd')
    r.postfit_phs = Numeric.zeros(r.numTOAs, 'd')
    r.postfit_sec = Numeric.zeros(r.numTOAs, 'd')
    r.orbit_phs = Numeric.zeros(r.numTOAs, 'd')
    r.bary_freq = Numeric.zeros(r.numTOAs, 'd')
    r.weight = Numeric.zeros(r.numTOAs, 'd')
    r.uncertainty = Numeric.zeros(r.numTOAs, 'd')
    r.prefit_phs = Numeric.zeros(r.numTOAs, 'd')
    for ii in range(r.numTOAs):
        struct = rf.fort_read("idddddddddi")
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
    if not Numeric.nonzero(r.orbit_phs): del r.orbit_phs
    if not Numeric.nonzero(r.bary_freq): del r.bary_freq
    if not Numeric.nonzero(r.weight): del r.weight
    r.prefit_sec = r.postfit_sec/r.postfit_phs*r.prefit_phs
    r.uncertainty *= 1.e-6 # Convert uncertainties in usec to sec
    return r

def read_residuals_64bit():
    r = residuals()
    rf = open("resid2.tmp")
    # Fortran format: 2 longs (Fortran crap) + 9 doubles = 80 bytes
    rf.seek(0, 2)
    r.numTOAs = rf.tell()/88
    rf.seek(0, 0)
    r.bary_TOA = Numeric.zeros(r.numTOAs, 'd')
    r.postfit_phs = Numeric.zeros(r.numTOAs, 'd')
    r.postfit_sec = Numeric.zeros(r.numTOAs, 'd')
    r.orbit_phs = Numeric.zeros(r.numTOAs, 'd')
    r.bary_freq = Numeric.zeros(r.numTOAs, 'd')
    r.weight = Numeric.zeros(r.numTOAs, 'd')
    r.uncertainty = Numeric.zeros(r.numTOAs, 'd')
    r.prefit_phs = Numeric.zeros(r.numTOAs, 'd')
    for ii in range(r.numTOAs):
	block = rf.read(88)
	struct = Numeric.fromstring(block, 'd')
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
    if not Numeric.nonzero(r.orbit_phs): del r.orbit_phs
    if not Numeric.nonzero(r.bary_freq): del r.bary_freq
    if not Numeric.nonzero(r.weight): del r.weight
    r.prefit_sec = r.postfit_sec/r.postfit_phs*r.prefit_phs
    r.uncertainty *= 1.e-6 # Convert uncertainties in usec to sec
    return r
