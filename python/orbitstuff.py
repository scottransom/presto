import random, presto, math, Numeric

def copyorb(old, new):
    """
    copyorb(old, new):
        Copy an orbitparams variable from 'old' to 'new'.
    """
    new.p = old.p
    new.e = old.e
    new.x = old.x
    new.w = old.w
    new.t = old.t
    new.wd = old.wd
    new.pd = old.pd
    return new

def fake_mspsr(companion='None', psrp=None, orbp=None,
               orbx=None, orbe=None, orbw=None, orbt=None):
    """
    fake_mspsr(companion='None'):
        Generate random pulsar parameters.
        Returns a psrparams structure.
        
        Keyword Arguments:
        companion -- the companion type ('WD', 'NS', 'BH', default = None)
            Note:  All of the following default to random values.
        psrp -- pulsar period in sec.
        orbp -- orbital period in sec.
        orbx -- projected orbital semi-major axis in lt-sec.
        orbe -- orbital eccentricity.
        orbw -- argument of periapsis in degrees.
        orbt -- time since last periapsis passage in sec.
    """
    global mpsr, mc, mwd, mns, mbh
    psr = presto.psrparams()
    psr.jname = 'Fake PSR'
    psr.bname = 'Fake PSR'
    psr.ntype = 0; psr.ra2000 = 0.0; psr.dec2000 = 0.0
    psr.dm = 0.0; psr.dist = 0.0
    psr.pd = 0.0; psr.pdd = 0.0
    psr.fd = 0.0; psr.fdd = 0.0
    if not psrp:  psr.p = random.uniform(0.002, 0.07)
    else: psr.p = psrp
    if psr.p < 0.03: psr.ntype = psr.ntype | 16
    psr.f = 1.0 / psr.p
    psr.fwhm = random.uniform(0.05, 0.5) * psr.p * 1000.0
    psr.timepoch = 0.0
    if companion=='WD':
        mwd = random.gauss(0.25, 0.10)
        mc = mwd;
        if not orbe: psr.orb.e = random.expovariate(80.0)
        else: psr.orb.e = orbe
        ctype = 'WD'
        psr.ntype = psr.ntype | 8
    elif companion=='NS':
        mns = random.gauss(1.35, 0.04)
        mc = mns;
        if not orbe: psr.orb.e = abs(random.uniform(0.1, 0.8))
        else: psr.orb.e = orbe
        ctype = 'NS'
        psr.ntype = psr.ntype | 8
    elif companion=='BH':
        mbh = random.gauss(6.0, 1.0)
        mc = mbh;
        if not orbe: psr.orb.e = abs(random.uniform(0.1, 0.8))
        else: psr.orb.e = orbe
        ctype = 'BH'
        psr.ntype = psr.ntype | 8
    else:
        psr.orb.p = 0.0; psr.orb.x = 0.0; psr.orb.e = 0.0;
        psr.orb.w = 0.0; psr.orb.t = 0.0; psr.orb.pd = 0.0;
        psr.orb.wd = 0.0;
        mc = 0.0; ctype = 'None'
    if ctype!='None':
        # The following is from Thorsett and Chakrabarty, 1988.
        mpsr = random.gauss(1.35, 0.04)
        inc = math.acos(random.uniform(0.0, 1.0))
        mf = (mc * math.sin(inc))**3 / (mc + mpsr)**2
        if not orbp: psr.orb.p = random.uniform(20.0, 600.0) * 60.0
        else: psr.orb.p = orbp
        if not orbx: psr.orb.x = (mf * psr.orb.p**2
                                  * 1.24764143192e-07)**(1.0 / 3.0)
        else: psr.orb.x = orbx
        if not orbw: psr.orb.w = random.uniform(0.0, 360.0)
        else: psr.orb.w = orbw
        if not orbt: psr.orb.t = random.uniform(0.0, psr.orb.p)
        else: psr.orb.t = orbt
    return psr

