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

import mIO, Numeric

class residuals:
    pass

def read_residuals():
    r = residuals()
    rf = mIO.binary_file("resid2.tmp")
    # Fortran format: 2 ints (Fortran crap) + 9 doubles = 80 bytes
    r.numTOAs = rf.size()/80
    r.bary_TOA = Numeric.zeros(r.numTOAs, 'd')
    r.postfit_phs = Numeric.zeros(r.numTOAs, 'd')
    r.postfit_sec = Numeric.zeros(r.numTOAs, 'd')
    r.orbit_phs = Numeric.zeros(r.numTOAs, 'd')
    r.bary_freq = Numeric.zeros(r.numTOAs, 'd')
    r.weight = Numeric.zeros(r.numTOAs, 'd')
    r.uncertainty = Numeric.zeros(r.numTOAs, 'd')
    r.prefit_phs = Numeric.zeros(r.numTOAs, 'd')
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
    return r
