from barycenter import *
import math
import fpformat

def posvec(ra, dec):
    """
    posvec(ra,dec):
        Convert an onbects RA and DEC into a unit vector pointing
        towards the source in B1950.0 coords.
        Arguments:
            'ra' is the RA of the object in B1950.0 coords (radians)
            'dec' is the DEC of the object in B1950.0 coords (radians)
    """
    return [math.cos(dec) * math.cos(ra),
            math.cos(dec) * math.sin(ra),
            math.sin(dec)]

def barwiz(filenm, dt_topo, numtimes, start_mjd, obs, ra, dec, freq, dm):
    """
    barwiz(filenm, dt_topo, numtimes, start_mjd, obs, ra, dec, freq, dm):
        Generate a file containing topocentric and barycentric time
        pairs.  The topocentric times should be UTC.  The barycentric
        times will be TDB.
        Arguments:
            'filenm' is the output filename
            'dt_topo' is the time between input samples
            'numtimes' is the number of conversions to perform
            'start_mjd' is the MJD (UTC) if the first topocentric time
            'obs' is a two character observatory code from TEMPO's obsys.dat
            'ra' is the J2000 right ascension of the observation
            'dec' is the J2000 declination of the observation
            'freq' is the observing frequency in hertz
            'dm' is the dispersion measure in pc/cm^3
    """
    bval = dm/2.410e-16
    site = obs_coords(obs)
    [r, d, dr, dd, p, v] = precess_J2000_to_B1950(ra, dec, 0.0, 0.0, 0.0, 0.0)
    pos = posvec(r,d)
    (fjdu, ijdu) = math.modf(start_mjd)
    ijdu = int(ijdu + 2400000)
    print ''
    print 'Writing Fortran barycenter timing code.'
    print ''
    f = open('baryrun.f', 'w')
    f.write('C|BARWIZ\n')
    f.write('       IMPLICIT REAL*8(A-H,O-Z)\n')
    f.write('       COMMON/CONST/ PI,TWOPI,SECDAY,CONVD,CONVS,AULTSC,VELC,EMRAT,OBLQ,\n')
    f.write('     1 GAUSS,RSCHW,AULTVL\n')
    f.write('       COMMON/CRDS/ RCB(6),RBE(6),RCE(6),RSE(6),RCA(6),RSA(6),RBA(6),\n')
    f.write('     1 S(6),PSID,EPSD,PC,PS,TSID,PRN(3,3),ATUT,UT1UT,ETAT,DTGR,TDIS,BCLT\n')
    f.write('       COMMON/INODTA/ IN,IOUT\n')
    f.write('       COMMON/OBSP/ SITE(3),POS(3),FREQ,BVAL\n')
    f.write('       DATA SITE/'+fpformat.fix(site[0],13)+','+fpformat.fix(site[1],13)+','+fpformat.fix(site[2],13)+'/\n')
    f.write('       DATA POS/'+fpformat.fix(pos[0],13)+','+fpformat.fix(pos[1],13)+','+fpformat.fix(pos[2],13)+'/\n')
    f.write('       DATA JDU,FRU /'+`ijdu`+','+fpformat.fix(fjdu,14)+'/\n')
    f.write("       OPEN(1,FILE='nbody740.bcd68t99',STATUS='UNKNOWN')\n")
    f.write("       OPEN(42,FILE='earthrot.dat',STATUS='UNKNOWN')\n")
    f.write("       OPEN(20,FILE='"+filenm+".tim',STATUS='UNKNOWN')\n")
    f.write('       CALL BARTIM(0,0D0,0,0D0)\n')
    f.write('       FREQ='+`freq`+'\n')
    f.write('       BVAL='+`bval`+'\n')
    f.write('       DFRAC=('+`dt_topo`+')/86400D0\n')
    f.write('       DO 10, CT=1,'+`numtimes`+'\n')
    f.write('          CALL BARTIM(JDU,FRU,JDC,FRC)\n')
    f.write('          WRITE(20,20) JDU-2400000,FRU,JDC-2400000,FRC\n')
    f.write(" 20       FORMAT(I5,F14.13,'  ',I5,F14.13)\n")
    f.write('          FRU=FRU+DFRAC\n')
    f.write('          IF (FRU .GE. 1.0D0) THEN\n')
    f.write('            FRU=FRU-1.0D0\n')
    f.write('            JDU=JDU+1\n')
    f.write('          ENDIF\n')
    f.write(' 10    CONTINUE\n')
    f.write('       CLOSE(1)\n')
    f.write('       CLOSE(20)\n')
    f.write('       CLOSE(42)\n')
    f.write('       STOP\n')
    f.write('       END\n')
    f.close()
#    print 'Compiling.'
#    print ''
#    cmd = 'f77 baryrun.f Barycenter/bartim2.f -o baryrun.exe'
#    spawn, cmd
#    print ''
#    print 'Running.'
#    print ''
#    spawn, 'baryrun.exe'
#    print ''
#    print 'Wrote JDU(I7).FRU(F14.13)" "JDC(I7).FRC(F14.13)'
#    print '   as formatted double precision in the file called:'
#    print '   '+filenm+'.tim'
#    print ''






