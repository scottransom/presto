import mIO, Numeric

def read_residuals():
    rf = mIO.binary_file("resid2.tmp")
    # Fortran format: 2 ints (Fortran crap) + 9 doubles = 80 bytes
    numTOAs = rf.size()/80
    pre = Numeric.zeros(numTOAs, 'd')
    post_sec = Numeric.zeros(numTOAs, 'd')
    post_phs = Numeric.zeros(numTOAs, 'd')
    weights = Numeric.zeros(numTOAs, 'd')
    orb_phs = Numeric.zeros(numTOAs, 'd')
    for ii in range(numTOAs):
        struct =rf.fort_read("ddddddddd")
        (pre[ii], post_sec[ii], post_phs[ii], weights[ii], orb_phs[ii]) = \
                  (struct[7], struct[2], struct[1], struct[5], struct[3])
    rf.close()
    return (pre, post_sec, post_phs, weights, orb_phs)
