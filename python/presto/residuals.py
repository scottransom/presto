from __future__ import print_function
from builtins import range
from builtins import object
import struct
import numpy as Num

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
#


class residuals(object):
    pass

def read_residuals(filename="resid2.tmp"):
    """
    read_residuals(filename="resid2.tmp"):
        Read a TEMPO1 style binary residuals file and return all the elements
            in a residuals 'class'.  The class instance will have an attribute
            called .numTOAs with the number of TOAs and up to 8 arrays with
            the following (as appropriate):
            .bary_TOA     Barycentric TOA (MJD)
            .uncertainty  TOA uncertainty (seconds)
            .bary_freq    Observing frequency (in barycenter frame)
            .prefit_phs   Prefit residual (pulse phase, from 0 to 1)
            .prefit_sec   Prefit residual (seconds)
            .postfit_phs  Postfit residual (pulse phase, from 0 to 1)
            .postfit_sec  Postfit residual (seconds)
            .orbit_phs    Orbital phase (where applicable)
            .weight       Weight of point in the fit
    """
    r = residuals()
    infile = open(filename, "rb")
    swapchar = '<' # this is little-endian (default)
    data = infile.read(8)
    test_int32 = struct.unpack(swapchar+"i", data[:4])[0]
    test_int64 = struct.unpack(swapchar+"q", data)[0]
    if ((test_int32 > 100 or test_int32 < 0) and
        (test_int64 > 100 or test_int64 < 0)):
        swapchar = '>' # this is big-endian
    if (test_int32 < 100 and test_int32 > 0):
        marktype = 'i'  # 32-bit int
        reclen = test_int32 + 2 * 4
    else:
        marktype = 'q'  # long long
        reclen = test_int64 + 2 * 8
    rectype = swapchar+marktype+9*'d'+marktype
    # print test_int32, test_int64, marktype, reclen, rectype
    infile.seek(0, 2) # position at file end
    filelen = infile.tell()
    if (filelen % reclen or
        not (reclen==struct.calcsize(rectype))):
        print("Warning:  possibly reading residuals incorrectly... don't understand record size")
    infile.seek(0, 0) # position at file start
    r.numTOAs = filelen // reclen
    r.bary_TOA = Num.zeros(r.numTOAs, 'd')
    r.postfit_phs = Num.zeros(r.numTOAs, 'd')
    r.postfit_sec = Num.zeros(r.numTOAs, 'd')
    r.orbit_phs = Num.zeros(r.numTOAs, 'd')
    r.bary_freq = Num.zeros(r.numTOAs, 'd')
    r.weight = Num.zeros(r.numTOAs, 'd')
    r.uncertainty = Num.zeros(r.numTOAs, 'd')
    r.prefit_phs = Num.zeros(r.numTOAs, 'd')
    for ii in range(r.numTOAs):
        rec = struct.unpack(rectype, infile.read(reclen))
        (r.bary_TOA[ii], 
         r.postfit_phs[ii],
         r.postfit_sec[ii], 
         r.orbit_phs[ii], 
         r.bary_freq[ii], 
         r.weight[ii], 
         r.uncertainty[ii], 
         r.prefit_phs[ii]) = (rec[1], rec[2], rec[3], rec[4], \
                              rec[5], rec[6], rec[7], rec[8])
    infile.close()
    if not Num.nonzero(r.orbit_phs): del r.orbit_phs
    if not Num.nonzero(r.bary_freq): del r.bary_freq
    if not Num.nonzero(r.weight): del r.weight
    r.prefit_sec = r.postfit_sec/r.postfit_phs*r.prefit_phs
    r.uncertainty *= 1.e-6 # Convert uncertainties in usec to sec
    return r
