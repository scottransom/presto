from builtins import range
## Automatically adapted for numpy Apr 14, 2006 by convertcode.py

import math
import numpy as Numeric 

# Most of these routines were taken from TEMPO v11.004 (arrtim.f, setup.f)

def convert_angle(inval, flag=1):
    """
    convert_angle(inval, flag=1):
        Converts a coded double to an angle based on the optional 'flag'
        'flag' = 1:  Input form = ddmmss.ss, Output in radians (default)
        'flag' = 2:  Input form = ddmmss.ss, Output in frac of 2pi
        'flag' = 3:  Input form = hhmmss.ss, Output in radians
        'flag' = 4:  Input form = hhmmss.ss, Output in frac of 2pi
    """
    twopi = 6.28318530717958647692528676655900576
    (im, id) = math.modf(inval / 10000.0)          # Integer deg/hours
    (s, im) = math.modf((im * 10000.0) / 100.0)    # Integer minutes
    s = s * 100.0                                  # seconds
    ang = (id + (im + s / 60.0) / 60.0) / 360.0    # Angle in frac of 2pi
    if (flag > 2): ang = ang * 15.0
    if (flag == 1 or flag == 3): ang = ang * twopi
    return ang

def hms2hours(hms):
    """
    hms2hours(hms):
        Converts an angle 'hms' expressed as hhmmss.ss into
        fractional hours.
    """
    return convert_angle(hms, 4) * 24.0

def dms2deg(dms):
    """
    dms2deg(dms):
        Converts an angle 'dms' expressed as ddmmss.ss into
        fractional degrees.
    """
    return convert_angle(dms, 2) * 360.0

def rad2dms(rad):
    """
    rad2dms(rad):
        Convert 'rad' radians into dd:mm:ss.sss format.
    """
    deg = rad * 180.0 / math.pi
    (mm, dd) = math.modf(deg)
    (ss, mm) = math.modf(mm * 60.0)
    ss = ss * 60.0
    id = abs(int(dd))
    if (abs(dd) < 10.0):
        if (dd < 0): d = '-0'+repr(id)
        else: d = '0'+repr(id)
    else:
        if (dd < 0): d = '-'+repr(id)
        else: d = repr(id)
    im = abs(int(mm))
    if (abs(mm) < 10): m = '0'+repr(im)
    else: m = repr(im)
    if (abs(ss) < 10): s = '0'+repr(abs(ss))
    else: s = repr(abs(ss))
    return d+':'+m+':'+s

def rad2hms(rad):
    """
    rad2hms(rad):
        Convert 'rad' radians into hh:mm:ss.sss format.
    """
    hours = rad * 12.0 / math.pi
    (mm, hh) = math.modf(hours)
    (ss, mm) = math.modf(mm * 60.0)
    ss = ss * 60.0
    ih = abs(int(hh))
    if (abs(hh) < 10.0):
        if (hh < 0): h = '-0'+repr(ih)
        else: h = '0'+repr(ih)
    else:
        if (hh < 0): h = '-'+repr(ih)
        else: h = repr(ih)
    im = abs(int(mm))
    if (abs(mm) < 10): m = '0'+repr(im)
    else: m = repr(im)
    if (abs(ss) < 10): s = '0'+repr(abs(ss))
    else: s = repr(abs(ss))
    return h+':'+m+':'+s


def geodetic2geocentcyl(lat, lon, elev):
    """
    geodetic2geocentcyl(lat, lon, elev):
        Return a list containing the Geocentric Cylindrical coords.
            'lat' is Geodetic latitude in degrees (ddmmss.ss)
            'long' is Geodetic west longitude in degrees (ddmmss.ss)
            'elev' is meters above mean sea level
    """
    rad_e = 6378140.0
    velc = 299792458.0
    flat = 1.0/298.257
    ault = 499.004786
    lon = convert_angle(lon)
    lat = convert_angle(lat)
    aa_c = 1.0 / math.sqrt(1.0 + (-2.0 + flat) * \
                           flat * math.sin(lat) * math.sin(lat))
    aa_arcf = (rad_e * aa_c + elev) * math.cos(lat)
    aa_arsf = (rad_e * (1.0 - flat) * (1.0 - flat) * \
               aa_c + elev) * math.sin(lat)
    hlt = math.atan2(aa_arsf, aa_arcf)
    erad = math.sqrt(aa_arcf * aa_arcf + aa_arsf * aa_arsf)
    hrd = erad / (velc * ault)
    site = []  # site is a list containing [r, z/velc, longitude in rad]
    site.append(hrd * math.cos(hlt) * ault)
    site.append(site[0] * math.tan(hlt))
    site.append(lon)
    return site

def xyz2geocentcyl(x, y, z):
    """
    xyz2geocentcyl(x, y, z):
        Return a list containing the Geocentric Cylindrical coords.
            'x', 'y', and 'z' are referenced to the Geocenter in m.
    """
    velc = 299792458.0
    ault = 499.004786
    erad = math.sqrt(x * x + y * y + z * z)
    hlt = math.asin(z / erad)
    lon = math.atan2(-y, x)
    hrd = erad / (velc * ault)
    site = []  # site is a list containing [r, z/velc, longitude in rad]
    site.append(hrd * math.cos(hlt) * ault)
    site.append(site[0] * math.tan(hlt))
    site.append(lon)
    return site

def obs_coords(observ):
    """
    obs_coords(observ):
        Return a list containing the Geocentric Cylindrical Coords for
        an observatory found in the TEMPO 'obsys.dat'.
            'observ' is the two letter observatory code
                from 'obsys.dat' (i.e. 'PK' = Parkes)
    """
    obs = {}
    obs['GB'] = [382645.48, 795054.53, 798.5, '', 'GREEN BANK']
    obs['QU'] = [422333.2, 722040.4, 306.0, '', 'QUABBIN']
    obs['AO'] = [2390490.0, -5564764.0, 1994727.0, 'XYZ', 'ARECIBO XYZ (JPL)']
    obs['HO'] = [-424818.0, -1472621.0, 50.0, '', 'Hobart, Tasmania']
    obs['PR'] = [402047.7, 743853.85, 43.0, '', 'PRINCETON']
    obs['VL'] = [-1601192.0, -5041981.4, 3554871.4, 'XYZ', 'VLA XYZ']
    obs['PK'] = [-330000.04, -1481542.00, 392.0, '', 'PARKES']
    obs['JB'] = [3822252.643, -153995.683, 5086051.443, 'XYZ', 'JODRELL BANK']
    obs['G3'] = [382546.30, 795056.36, 893.7, '', 'GB 300FT']
    obs['G1'] = [382615.409, 795009.613, 880.87, '', 'GB 140FT']
    obs['G8'] = [382545.90, 795036.87, 835.8, '', 'GB 85-3']
    obs['V2'] = [340443.497, 1073703.819, 2124.0, '', 'VLA SITE']
    obs['BO'] = [443118.48, -113848.48, 25.0, '', 'NORTHERN CROSS']
    obs['MO'] = [-352219.00, -1492525.00, 500.0, '', 'MOST']
    obs['NC'] = [4324165.81, 165927.11, 4670132.83, 'XYZ', 'Nancay']
    obs['EF'] = [4033949.5, 486989.4, 4900430.8, 'XYZ', 'Effelsberg']
    obs['JB'] = [531412.0, 21824.0, 78.0, '', 'JODRELL BANK']
    obs['FB'] = [332235.0, 1171501.0, 0.0, '', 'Fallbrook']
    obs['MT'] = [314119.6, 1105304.4, 2606.0, '', 'MMT']
    if obs[observ][3] == 'XYZ':
        return xyz2geocentcyl(obs[observ][0], obs[observ][1], \
                              obs[observ][2])
    else:
        return geodetic2geocentcyl(obs[observ][0], obs[observ][1], \
                                   obs[observ][2])

def precess_J2000_to_B1950(ra, dec, rapm=0.0, decpm=0.0, par=0.0, vel=0.0): 
    """
    precess_J2000_to_B1950(ra, dec, rapm, decpm, par, vel): 
        Precess a set of J2000.0 FK5 coords to epoch B1950.0 FK4.
        Return a list containing the B1950 version of the arguments.
            'ra' is the J2000 right ascension (hhmmss.ssss)
            'dec' is the J2000 declination (ddmmss.ssss)
            'rapm' is the J2000 RA proper motion in rad/Julian Year (0.0)
            'decpm' is the J2000 DEC proper motion in rad/Julian Year (0.0)
            'par' is the parallax in arcsec (0.0)
            'vel' is the radial velocity in km/s (+ = away from us) (0.0)
        Note: Parenthesized values at the ends of the above lines are the
              default values.
    """
    # This is converted from the SLALIB routine fk524.f
    # This routine converts stars from the new, IAU 1976, FK5, Fricke
    # system, to the old, Bessel-Newcomb, FK4 system.  The precepts
    # of Smith et al (Ref 1) are followed, using the implementation
    # by Yallop et al (Ref 2) of a matrix method due to Standish.
    # Kinoshita's development of Andoyer's post-Newcomb precession is
    # used.  The numerical constants from Seidelmann et al (Ref 3) are
    # used canonically.
    twopi = 6.283185307179586476925287
    # Radians per year to arcsec per century 
    pmf = 100.0 * 60.0 * 60.0 * 360.0 / twopi
    tiny = 1.0e-30
    # Km per sec to AU per tropical century
    #    = 86400 * 36524.2198782 / 149597870
    vf = 21.095
    a = Numeric.array([-1.62557e-6, -0.31919e-6, -0.13843e-6, \
                       +1.245e-3, -1.580e-3, -0.659e-3])
    emi = Numeric.array([[+0.9999256795, -0.0111814828, -0.0048590040,
                          -0.000551, -0.238560, +0.435730],
                         [+0.0111814828, +0.9999374849, -0.0000271557,
                          +0.238509, -0.002667, -0.008541],
                         [+0.0048590039, -0.0000271771, +0.9999881946,
                          -0.435614, +0.012254, +0.002117],
                         [-2.42389840e-6, +2.710544e-8, +1.177742e-8,
                          +0.99990432, -0.01118145, -0.00485852],
                         [-2.710544e-8, -242392702e-6, +6.585e-11,
                          +0.01118145, +0.99991613, -0.00002716],
                         [-1.177742e-8, +6.585e-11, -2.42404995e-6,
                          +0.00485852, -0.00002717,+0.99996684]])
    r = convert_angle(ra, 3)
    d = convert_angle(dec, 1)
    ur = rapm * pmf
    ud = decpm * pmf
    px = par
    rv = vel
    sr = math.sin(r)
    cr = math.cos(r)
    sd = math.sin(d)
    cd = math.cos(d)
    w = vf * rv * px
    x = cr * cd
    y = sr * cd
    z = sd
    v1 = Numeric.array([x, y, z, -ur * y - cr * sd * ud + w * x,
                        ur * x - sr * sd * ud + w * y, cd * ud + w * z])
    # convert position+velocity vector to BN system
    v2 = Numeric.zeros(6, 'd')
    for i in range(6):
        w = 0.0
        for j in range(6):
            w = w + emi[j][i] * v1[j]
        v2[i] = w
    # position vector components and magnitude
    x = v2[0]
    y = v2[1]
    z = v2[2]
    rxyz = math.sqrt(x * x + y * y + z * z)
    # apply e-terms to position
    w = x * a[0] + y * a[1] + z * a[2]
    x = x + a[0] * rxyz - w * x
    y = y + a[1] * rxyz - w * y
    z = z + a[2] * rxyz - w * z
    # recompute magnitude
    rxyz = math.sqrt(x * x + y * y + z * z)
    # apply e-terms to both position and velocity
    x = v2[0]
    y = v2[1]
    z = v2[2]
    w = x * a[0] + y * a[1] + z * a[2]
    wd = x * a[3] + y * a[4] + z * a[5]
    x = x + a[0] * rxyz - w * x
    y = y + a[1] * rxyz - w * y
    z = z + a[2] * rxyz - w * z
    xd = v2[3] + a[3] * rxyz - wd * x
    yd = v2[4] + a[4] * rxyz - wd * y
    zd = v2[5] + a[5] * rxyz - wd * z
    # convert to spherical
    rxysq = x * x + y * y
    rxy = math.sqrt(rxysq)
    if (x==0.0 and y==0.0): r = 0.0
    else:
        r = math.atan2(y, x)
        if (r < 0.0): r = r + twopi
    d = math.atan2(z, rxy)
    if (rxy > tiny):
        ur = (x * yd - y * xd) / rxysq
        ud = (zd * rxysq - z * (x * xd + y * yd)) / ((rxysq + z * z) * rxy)
    # radial velocity and parallax
    if (px > tiny):
        rv = (x * xd + y * yd + z * zd) / (px * vf * rxyz)
        px = px / rxyz
    return [r, d, ur / pmf, ud / pmf, px, rv]


# Most of the following time formulae are from
#    http://legacy.gsfc.nasa.gov/docs/xte/abc/time_tutorial.html

def TAI_minus_UTC(mjd):
    """
    TAI_minus_UTC(mjd):
        Return the difference between TAI (International Atomic Time)
        and UTC in seconds at a specified MJD.
            'mjd' is the Modified Julian Date UTC
    """
    LeapSecMJD = Numeric.array([41499, 41683, 42048, 42413, 42778, \
                                43144, 43509, 43874, 44239, 44786, \
                                45151, 45516, 46247, 47161, 47892, \
                                48257, 48804, 49169, 49534, 50083, \
                                50630, 51179])
    TAI_UTC_diffs = Numeric.arange(len(LeapSecMJD)+1)+10
    TAI_minus_UTC = TAI_UTC_diffs[Numeric.searchsorted(LeapSecMJD, mjd)]
    return TAI_minus_UTC

def TT_minus_UTC(mjd):
    """
    TT_minus_UTC(mjd):
        Return the difference between TT (Terrestrial Dynamic Time)
        and UTC in seconds at a specified MJD.  TT used to be called
        ET (Ephemeris Time).
            'mjd' is the Modified Julian Date UTC
    """
    TT_minus_TAI = 32.184
    TT_minus_UTC = TT_minus_TAI + TAI_minus_UTC(mjd)
    return TT_minus_UTC

def TDB_minus_UTC(mjd):
    """
    TDB_minus_UTC(mjd):
        Return the difference between TDB (Barycentric Dynamic Time)
        and UTC in seconds at a specified MJD.
            'mjd' is the Modified Julian Date UTC
    """
    g = (357.53 + 0.9856003 * (mjd - 51544.5)) * 0.017453292519943295769
    TDB_minus_TT = 0.001658 * math.sin(g) + 0.000014 * math.sin(2.0 * g)
    TDB_minus_UTC = TDB_minus_TT + TT_minus_UTC(mjd)
    return TDB_minus_UTC
