      SUBROUTINE sla_PLANET (DATE, NP, PV, JSTAT)
*+
*     - - - - - - -
*      P L A N E T
*     - - - - - - -
*
*  Approximate heliocentric position and velocity of a specified
*  major planet.
*
*  Given:
*     DATE      d      Modified Julian Date (JD - 2400000.5)
*     NP        i      planet (1=Mercury, 2=Venus, 3=EMB ... 9=Pluto)
*
*  Returned:
*     PV        d(6)   heliocentric x,y,z,xdot,ydot,zdot, J2000
*                                           equatorial triad (AU,AU/s)
*     JSTAT     i      status: +1 = warning: date out of range
*                               0 = OK
*                              -1 = illegal NP (outside 1-9)
*                              -2 = solution didn't converge
*
*  Called:  sla_PLANEL
*
*  Notes
*
*  1  The epoch, DATE, is in the TDB timescale and is a Modified
*     Julian Date (JD-2400000.5).
*
*  2  The reference frame is equatorial and is with respect to the
*     mean equinox and ecliptic of epoch J2000.
*
*  3  If an NP value outside the range 1-9 is supplied, an error
*     status (JSTAT = -1) is returned and the PV vector set to zeroes.
*
*  4  The algorithm for obtaining the mean elements of the planets
*     from Mercury to Neptune is due to J.L. Simon, P. Bretagnon,
*     J. Chapront, M. Chapront-Touze, G. Francou and J. Laskar
*     (Bureau des Longitudes, Paris).  The (completely different)
*     algorithm for calculating the ecliptic coordinates of Pluto
*     is by Meeus.
*
*  5  Comparisons of the present routine with the JPL DE200 ephemeris
*     give the following RMS errors over the interval 1960-2025:
*
*                      position (km)     speed (metre/sec)
*
*        Mercury            334               0.437
*        Venus             1060               0.855
*        EMB               2010               0.815
*        Mars              7690               1.98
*        Jupiter          71700               7.70
*        Saturn          199000              19.4
*        Uranus          564000              16.4
*        Neptune         158000              14.4
*        Pluto            36400               0.137
*
*     From comparisons with DE102, Simon et al quote the following
*     longitude accuracies over the interval 1800-2200:
*
*        Mercury                 4"
*        Venus                   5"
*        EMB                     6"
*        Mars                   17"
*        Jupiter                71"
*        Saturn                 81"
*        Uranus                 86"
*        Neptune                11"
*
*     In the case of Pluto, Meeus quotes an accuracy of 0.6 arcsec
*     in longitude and 0.2 arcsec in latitude for the period
*     1885-2099.
*
*     For all except Pluto, over the period 1000-3000 the accuracy
*     is better than 1.5 times that over 1800-2200.  Outside the
*     period 1000-3000 the accuracy declines.  For Pluto the
*     accuracy declines rapidly outside the period 1885-2099.
*     Outside these ranges (1885-2099 for Pluto, 1000-3000 for
*     the rest) a "date out of range" warning status (JSTAT=+1)
*     is returned.
*
*  6  The algorithms for (i) Mercury through Neptune and (ii) Pluto
*     are completely independent.  In the Mercury through Neptune
*     case, the present SLALIB implementation differs from the
*     original Simon et al Fortran code in the following respects.
*
*     *  The date is supplied as a Modified Julian Date rather
*        than a Julian Date (MJD = JD - 2400000.5).
*
*     *  The result is returned only in equatorial Cartesian form;
*        the ecliptic longitude, latitude and radius vector are not
*        returned.
*
*     *  The velocity is in AU per second, not AU per day.
*
*     *  Different error/warning status values are used.
*
*     *  Kepler's equation is not solved inline.
*
*     *  Polynomials in T are nested to minimize rounding errors.
*
*     *  Explicit double-precision constants are used to avoid
*        mixed-mode expressions.
*
*     *  There are other, cosmetic, changes to comply with
*        Starlink/SLALIB style guidelines.
*
*     None of the above changes affects the result significantly.
*
*  7  For NP=3 the result is for the Earth-Moon Barycentre.  To
*     obtain the heliocentric position and velocity of the Earth,
*     either use the SLALIB routine sla_EVP (or sla_EPV) or call
*     sla_DMOON and subtract 0.012150581 times the geocentric Moon
*     vector from the EMB vector produced by the present routine.
*     (The Moon vector should be precessed to J2000 first, but this
*     can be omitted for modern epochs without introducing significant
*     inaccuracy.)
*
*  References:  Simon et al., Astron. Astrophys. 282, 663 (1994).
*               Meeus, Astronomical Algorithms, Willmann-Bell (1991).
*
*  This revision:  19 June 2004
*
*  Copyright (C) 2004 P.T.Wallace.  All rights reserved.
*
*  License:
*    This program is free software; you can redistribute it and/or modify
*    it under the terms of the GNU General Public License as published by
*    the Free Software Foundation; either version 2 of the License, or
*    (at your option) any later version.
*
*    This program is distributed in the hope that it will be useful,
*    but WITHOUT ANY WARRANTY; without even the implied warranty of
*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*    GNU General Public License for more details.
*
*    You should have received a copy of the GNU General Public License
*    along with this program (see SLA_CONDITIONS); if not, write to the 
*    Free Software Foundation, Inc., 59 Temple Place, Suite 330, 
*    Boston, MA  02111-1307  USA
*
*-

      IMPLICIT NONE

      DOUBLE PRECISION DATE
      INTEGER NP
      DOUBLE PRECISION PV(6)
      INTEGER JSTAT

*  2Pi, deg to radians, arcsec to radians
      DOUBLE PRECISION D2PI,D2R,AS2R
      PARAMETER (D2PI=6.283185307179586476925286766559D0,
     :           D2R=0.017453292519943295769236907684886D0,
     :           AS2R=4.848136811095359935899141023579D-6)

*  Gaussian gravitational constant (exact)
      DOUBLE PRECISION GCON
      PARAMETER (GCON=0.01720209895D0)

*  Seconds per Julian century
      DOUBLE PRECISION SPC
      PARAMETER (SPC=36525D0*86400D0)

*  Sin and cos of J2000 mean obliquity (IAU 1976)
      DOUBLE PRECISION SE,CE
      PARAMETER (SE=0.3977771559319137D0,
     :           CE=0.9174820620691818D0)

      INTEGER I,J,IJSP(3,43)
      DOUBLE PRECISION AMAS(8),A(3,8),DLM(3,8),E(3,8),
     :                 PI(3,8),DINC(3,8),OMEGA(3,8),
     :                 DKP(9,8),CA(9,8),SA(9,8),
     :                 DKQ(10,8),CLO(10,8),SLO(10,8),
     :                 T,DA,DE,DPE,DI,DO,DMU,ARGA,ARGL,DM,
     :                 AB(2,3,43),DJ0,DS0,DP0,DL0,DLD0,DB0,DR0,
     :                 DJ,DS,DP,DJD,DSD,DPD,WLBR(3),WLBRD(3),
     :                 WJ,WS,WP,AL,ALD,SAL,CAL,
     :                 AC,BC,DL,DLD,DB,DBD,DR,DRD,
     :                 SL,CL,SB,CB,SLCB,CLCB,X,Y,Z,XD,YD,ZD

*  -----------------------
*  Mercury through Neptune
*  -----------------------

*  Planetary inverse masses
      DATA AMAS / 6023600D0,408523.5D0,328900.5D0,3098710D0,
     :            1047.355D0,3498.5D0,22869D0,19314D0 /

*
*  Tables giving the mean Keplerian elements, limited to T**2 terms:
*
*         A       semi-major axis (AU)
*         DLM     mean longitude (degree and arcsecond)
*         E       eccentricity
*         PI      longitude of the perihelion (degree and arcsecond)
*         DINC    inclination (degree and arcsecond)
*         OMEGA   longitude of the ascending node (degree and arcsecond)
*
      DATA A /
     :  0.3870983098D0,             0D0,      0D0,
     :  0.7233298200D0,             0D0,      0D0,
     :  1.0000010178D0,             0D0,      0D0,
     :  1.5236793419D0,           3D-10,      0D0,
     :  5.2026032092D0,       19132D-10,  -39D-10,
     :  9.5549091915D0, -0.0000213896D0,  444D-10,
     : 19.2184460618D0,       -3716D-10,  979D-10,
     : 30.1103868694D0,      -16635D-10,  686D-10 /
*
      DATA DLM /
     : 252.25090552D0, 5381016286.88982D0,  -1.92789D0,
     : 181.97980085D0, 2106641364.33548D0,   0.59381D0,
     : 100.46645683D0, 1295977422.83429D0,  -2.04411D0,
     : 355.43299958D0,  689050774.93988D0,   0.94264D0,
     :  34.35151874D0,  109256603.77991D0, -30.60378D0,
     :  50.07744430D0,   43996098.55732D0,  75.61614D0,
     : 314.05500511D0,   15424811.93933D0,  -1.75083D0,
     : 304.34866548D0,    7865503.20744D0,   0.21103D0/
*
      DATA E /
     : 0.2056317526D0,  0.0002040653D0,      -28349D-10,
     : 0.0067719164D0, -0.0004776521D0,       98127D-10,
     : 0.0167086342D0, -0.0004203654D0, -0.0000126734D0,
     : 0.0934006477D0,  0.0009048438D0,      -80641D-10,
     : 0.0484979255D0,  0.0016322542D0, -0.0000471366D0,
     : 0.0555481426D0, -0.0034664062D0, -0.0000643639D0,
     : 0.0463812221D0, -0.0002729293D0,  0.0000078913D0,
     : 0.0094557470D0,  0.0000603263D0,            0D0 /
*
      DATA PI /
     :  77.45611904D0,  5719.11590D0,   -4.83016D0,
     : 131.56370300D0,   175.48640D0, -498.48184D0,
     : 102.93734808D0, 11612.35290D0,   53.27577D0,
     : 336.06023395D0, 15980.45908D0,  -62.32800D0,
     :  14.33120687D0,  7758.75163D0,  259.95938D0,
     :  93.05723748D0, 20395.49439D0,  190.25952D0,
     : 173.00529106D0,  3215.56238D0,  -34.09288D0,
     :  48.12027554D0,  1050.71912D0,   27.39717D0 /
*
      DATA DINC /
     : 7.00498625D0, -214.25629D0,   0.28977D0,
     : 3.39466189D0,  -30.84437D0, -11.67836D0,
     :          0D0,  469.97289D0,  -3.35053D0,
     : 1.84972648D0, -293.31722D0,  -8.11830D0,
     : 1.30326698D0,  -71.55890D0,  11.95297D0,
     : 2.48887878D0,   91.85195D0, -17.66225D0,
     : 0.77319689D0,  -60.72723D0,   1.25759D0,
     : 1.76995259D0,    8.12333D0,   0.08135D0 /
*
      DATA OMEGA /
     :  48.33089304D0,  -4515.21727D0,  -31.79892D0,
     :  76.67992019D0, -10008.48154D0,  -51.32614D0,
     : 174.87317577D0,  -8679.27034D0,   15.34191D0,
     :  49.55809321D0, -10620.90088D0, -230.57416D0,
     : 100.46440702D0,   6362.03561D0,  326.52178D0,
     : 113.66550252D0,  -9240.19942D0,  -66.23743D0,
     :  74.00595701D0,   2669.15033D0,  145.93964D0,
     : 131.78405702D0,   -221.94322D0,   -0.78728D0 /
*
*  Tables for trigonometric terms to be added to the mean elements
*  of the semi-major axes.
*
      DATA DKP /
     : 69613, 75645, 88306, 59899, 15746, 71087, 142173,  3086,    0,
     : 21863, 32794, 26934, 10931, 26250, 43725,  53867, 28939,    0,
     : 16002, 21863, 32004, 10931, 14529, 16368,  15318, 32794,    0,
     : 6345,   7818, 15636,  7077,  8184, 14163,   1107,  4872,    0,
     : 1760,   1454,  1167,   880,   287,  2640,     19,  2047, 1454,
     :  574,      0,   880,   287,    19,  1760,   1167,   306,  574,
     :  204,      0,   177,  1265,     4,   385,    200,   208,  204,
     :    0,    102,   106,     4,    98,  1367,    487,   204,    0 /
*
      DATA CA /
     :      4,    -13,    11,    -9,    -9,    -3,    -1,     4,    0,
     :   -156,     59,   -42,     6,    19,   -20,   -10,   -12,    0,
     :     64,   -152,    62,    -8,    32,   -41,    19,   -11,    0,
     :    124,    621,  -145,   208,    54,   -57,    30,    15,    0,
     : -23437,  -2634,  6601,  6259, -1507, -1821,  2620, -2115,-1489,
     :  62911,-119919, 79336, 17814,-24241, 12068,  8306, -4893, 8902,
     : 389061,-262125,-44088,  8387,-22976, -2093,  -615, -9720, 6633,
     :-412235,-157046,-31430, 37817, -9740,   -13, -7449,  9644,    0 /
*
      DATA SA /
     :     -29,    -1,     9,     6,    -6,     5,     4,     0,    0,
     :     -48,  -125,   -26,   -37,    18,   -13,   -20,    -2,    0,
     :    -150,   -46,    68,    54,    14,    24,   -28,    22,    0,
     :    -621,   532,  -694,   -20,   192,   -94,    71,   -73,    0,
     :  -14614,-19828, -5869,  1881, -4372, -2255,   782,   930,  913,
     :  139737,     0, 24667, 51123, -5102,  7429, -4095, -1976,-9566,
     : -138081,     0, 37205,-49039,-41901,-33872,-27037,-12474,18797,
     :       0, 28492,133236, 69654, 52322,-49577,-26430, -3593,    0 /
*
*  Tables giving the trigonometric terms to be added to the mean
*  elements of the mean longitudes.
*
      DATA DKQ /
     :  3086, 15746, 69613, 59899, 75645, 88306, 12661, 2658,  0,   0,
     : 21863, 32794, 10931,    73,  4387, 26934,  1473, 2157,  0,   0,
     :    10, 16002, 21863, 10931,  1473, 32004,  4387,   73,  0,   0,
     :    10,  6345,  7818,  1107, 15636,  7077,  8184,  532, 10,   0,
     :    19,  1760,  1454,   287,  1167,   880,   574, 2640, 19,1454,
     :    19,   574,   287,   306,  1760,    12,    31,   38, 19, 574,
     :     4,   204,   177,     8,    31,   200,  1265,  102,  4, 204,
     :     4,   102,   106,     8,    98,  1367,   487,  204,  4, 102 /
*
      DATA CLO /
     :     21,   -95, -157,   41,   -5,   42,   23,   30,     0,    0,
     :   -160,  -313, -235,   60,  -74,  -76,  -27,   34,     0,    0,
     :   -325,  -322,  -79,  232,  -52,   97,   55,  -41,     0,    0,
     :   2268,  -979,  802,  602, -668,  -33,  345,  201,   -55,    0,
     :   7610, -4997,-7689,-5841,-2617, 1115, -748, -607,  6074,  354,
     : -18549, 30125,20012, -730,  824,   23, 1289, -352,-14767,-2062,
     :-135245,-14594, 4197,-4030,-5630,-2898, 2540, -306,  2939, 1986,
     :  89948,  2103, 8963, 2695, 3682, 1648,  866, -154, -1963, -283 /
*
      DATA SLO /
     :   -342,   136,  -23,   62,   66,  -52,  -33,   17,     0,    0,
     :    524,  -149,  -35,  117,  151,  122,  -71,  -62,     0,    0,
     :   -105,  -137,  258,   35, -116,  -88, -112,  -80,     0,    0,
     :    854,  -205, -936, -240,  140, -341,  -97, -232,   536,    0,
     : -56980,  8016, 1012, 1448,-3024,-3710,  318,  503,  3767,  577,
     : 138606,-13478,-4964, 1441,-1319,-1482,  427, 1236, -9167,-1918,
     :  71234,-41116, 5334,-4935,-1848,   66,  434,-1748,  3780, -701,
     : -47645, 11647, 2166, 3194,  679,    0, -244, -419, -2531,   48 /

*  -----
*  Pluto
*  -----

*
*  Coefficients for fundamental arguments:  mean longitudes
*  (degrees) and mean rate of change of longitude (degrees per
*  Julian century) for Jupiter, Saturn and Pluto
*
      DATA DJ0, DJD / 34.35D0, 3034.9057D0 /
      DATA DS0, DSD / 50.08D0, 1222.1138D0 /
      DATA DP0, DPD / 238.96D0, 144.9600D0 /

*  Coefficients for latitude, longitude, radius vector
      DATA DL0,DLD0 / 238.956785D0, 144.96D0 /
      DATA DB0 / -3.908202D0 /
      DATA DR0 / 40.7247248D0 /

*
*  Coefficients for periodic terms (Meeus's Table 36.A)
*
*  The coefficients for term n in the series are:
*
*    IJSP(1,n)     J
*    IJSP(2,n)     S
*    IJSP(3,n)     P
*    AB(1,1,n)     longitude sine (degrees)
*    AB(2,1,n)     longitude cosine (degrees)
*    AB(1,2,n)     latitude sine (degrees)
*    AB(2,2,n)     latitude cosine (degrees)
*    AB(1,3,n)     radius vector sine (AU)
*    AB(2,3,n)     radius vector cosine (AU)
*
      DATA (IJSP(I, 1),I=1,3),((AB(J,I, 1),J=1,2),I=1,3) /
     :                             0,  0,  1,
     :            -19798886D-6,  19848454D-6,
     :             -5453098D-6, -14974876D-6,
     :             66867334D-7,  68955876D-7 /
      DATA (IJSP(I, 2),I=1,3),((AB(J,I, 2),J=1,2),I=1,3) /
     :                             0,  0,  2,
     :               897499D-6,  -4955707D-6,
     :              3527363D-6,   1672673D-6,
     :            -11826086D-7,   -333765D-7 /
      DATA (IJSP(I, 3),I=1,3),((AB(J,I, 3),J=1,2),I=1,3) /
     :                             0,  0,  3,
     :               610820D-6,   1210521D-6,
     :             -1050939D-6,    327763D-6,
     :              1593657D-7,  -1439953D-7 /
      DATA (IJSP(I, 4),I=1,3),((AB(J,I, 4),J=1,2),I=1,3) /
     :                             0,  0,  4,
     :              -341639D-6,   -189719D-6,
     :               178691D-6,   -291925D-6,
     :               -18948D-7,    482443D-7 /
      DATA (IJSP(I, 5),I=1,3),((AB(J,I, 5),J=1,2),I=1,3) /
     :                             0,  0,  5,
     :               129027D-6,    -34863D-6,
     :                18763D-6,    100448D-6,
     :               -66634D-7,    -85576D-7 /
      DATA (IJSP(I, 6),I=1,3),((AB(J,I, 6),J=1,2),I=1,3) /
     :                             0,  0,  6,
     :               -38215D-6,     31061D-6,
     :               -30594D-6,    -25838D-6,
     :                30841D-7,     -5765D-7 /
      DATA (IJSP(I, 7),I=1,3),((AB(J,I, 7),J=1,2),I=1,3) /
     :                             0,  1, -1,
     :                20349D-6,     -9886D-6,
     :                 4965D-6,     11263D-6,
     :                -6140D-7,     22254D-7 /
      DATA (IJSP(I, 8),I=1,3),((AB(J,I, 8),J=1,2),I=1,3) /
     :                             0,  1,  0,
     :                -4045D-6,     -4904D-6,
     :                  310D-6,      -132D-6,
     :                 4434D-7,      4443D-7 /
      DATA (IJSP(I, 9),I=1,3),((AB(J,I, 9),J=1,2),I=1,3) /
     :                             0,  1,  1,
     :                -5885D-6,     -3238D-6,
     :                 2036D-6,      -947D-6,
     :                -1518D-7,       641D-7 /
      DATA (IJSP(I,10),I=1,3),((AB(J,I,10),J=1,2),I=1,3) /
     :                             0,  1,  2,
     :                -3812D-6,      3011D-6,
     :                   -2D-6,      -674D-6,
     :                   -5D-7,       792D-7 /
      DATA (IJSP(I,11),I=1,3),((AB(J,I,11),J=1,2),I=1,3) /
     :                             0,  1,  3,
     :                 -601D-6,      3468D-6,
     :                 -329D-6,      -563D-6,
     :                  518D-7,       518D-7 /
      DATA (IJSP(I,12),I=1,3),((AB(J,I,12),J=1,2),I=1,3) /
     :                             0,  2, -2,
     :                 1237D-6,       463D-6,
     :                  -64D-6,        39D-6,
     :                  -13D-7,      -221D-7 /
      DATA (IJSP(I,13),I=1,3),((AB(J,I,13),J=1,2),I=1,3) /
     :                             0,  2, -1,
     :                 1086D-6,      -911D-6,
     :                  -94D-6,       210D-6,
     :                  837D-7,      -494D-7 /
      DATA (IJSP(I,14),I=1,3),((AB(J,I,14),J=1,2),I=1,3) /
     :                             0,  2,  0,
     :                  595D-6,     -1229D-6,
     :                   -8D-6,      -160D-6,
     :                 -281D-7,       616D-7 /
      DATA (IJSP(I,15),I=1,3),((AB(J,I,15),J=1,2),I=1,3) /
     :                             1, -1,  0,
     :                 2484D-6,      -485D-6,
     :                 -177D-6,       259D-6,
     :                  260D-7,      -395D-7 /
      DATA (IJSP(I,16),I=1,3),((AB(J,I,16),J=1,2),I=1,3) /
     :                             1, -1,  1,
     :                  839D-6,     -1414D-6,
     :                   17D-6,       234D-6,
     :                 -191D-7,      -396D-7 /
      DATA (IJSP(I,17),I=1,3),((AB(J,I,17),J=1,2),I=1,3) /
     :                             1,  0, -3,
     :                 -964D-6,      1059D-6,
     :                  582D-6,      -285D-6,
     :                -3218D-7,       370D-7 /
      DATA (IJSP(I,18),I=1,3),((AB(J,I,18),J=1,2),I=1,3) /
     :                             1,  0, -2,
     :                -2303D-6,     -1038D-6,
     :                 -298D-6,       692D-6,
     :                 8019D-7,     -7869D-7 /
      DATA (IJSP(I,19),I=1,3),((AB(J,I,19),J=1,2),I=1,3) /
     :                             1,  0, -1,
     :                 7049D-6,       747D-6,
     :                  157D-6,       201D-6,
     :                  105D-7,     45637D-7 /
      DATA (IJSP(I,20),I=1,3),((AB(J,I,20),J=1,2),I=1,3) /
     :                             1,  0,  0,
     :                 1179D-6,      -358D-6,
     :                  304D-6,       825D-6,
     :                 8623D-7,      8444D-7 /
      DATA (IJSP(I,21),I=1,3),((AB(J,I,21),J=1,2),I=1,3) /
     :                             1,  0,  1,
     :                  393D-6,       -63D-6,
     :                 -124D-6,       -29D-6,
     :                 -896D-7,      -801D-7 /
      DATA (IJSP(I,22),I=1,3),((AB(J,I,22),J=1,2),I=1,3) /
     :                             1,  0,  2,
     :                  111D-6,      -268D-6,
     :                   15D-6,         8D-6,
     :                  208D-7,      -122D-7 /
      DATA (IJSP(I,23),I=1,3),((AB(J,I,23),J=1,2),I=1,3) /
     :                             1,  0,  3,
     :                  -52D-6,      -154D-6,
     :                    7D-6,        15D-6,
     :                 -133D-7,        65D-7 /
      DATA (IJSP(I,24),I=1,3),((AB(J,I,24),J=1,2),I=1,3) /
     :                             1,  0,  4,
     :                  -78D-6,       -30D-6,
     :                    2D-6,         2D-6,
     :                  -16D-7,         1D-7 /
      DATA (IJSP(I,25),I=1,3),((AB(J,I,25),J=1,2),I=1,3) /
     :                             1,  1, -3,
     :                  -34D-6,       -26D-6,
     :                    4D-6,         2D-6,
     :                  -22D-7,         7D-7 /
      DATA (IJSP(I,26),I=1,3),((AB(J,I,26),J=1,2),I=1,3) /
     :                             1,  1, -2,
     :                  -43D-6,         1D-6,
     :                    3D-6,         0D-6,
     :                   -8D-7,        16D-7 /
      DATA (IJSP(I,27),I=1,3),((AB(J,I,27),J=1,2),I=1,3) /
     :                             1,  1, -1,
     :                  -15D-6,        21D-6,
     :                    1D-6,        -1D-6,
     :                    2D-7,         9D-7 /
      DATA (IJSP(I,28),I=1,3),((AB(J,I,28),J=1,2),I=1,3) /
     :                             1,  1,  0,
     :                   -1D-6,        15D-6,
     :                    0D-6,        -2D-6,
     :                   12D-7,         5D-7 /
      DATA (IJSP(I,29),I=1,3),((AB(J,I,29),J=1,2),I=1,3) /
     :                             1,  1,  1,
     :                    4D-6,         7D-6,
     :                    1D-6,         0D-6,
     :                    1D-7,        -3D-7 /
      DATA (IJSP(I,30),I=1,3),((AB(J,I,30),J=1,2),I=1,3) /
     :                             1,  1,  3,
     :                    1D-6,         5D-6,
     :                    1D-6,        -1D-6,
     :                    1D-7,         0D-7 /
      DATA (IJSP(I,31),I=1,3),((AB(J,I,31),J=1,2),I=1,3) /
     :                             2,  0, -6,
     :                    8D-6,         3D-6,
     :                   -2D-6,        -3D-6,
     :                    9D-7,         5D-7 /
      DATA (IJSP(I,32),I=1,3),((AB(J,I,32),J=1,2),I=1,3) /
     :                             2,  0, -5,
     :                   -3D-6,         6D-6,
     :                    1D-6,         2D-6,
     :                    2D-7,        -1D-7 /
      DATA (IJSP(I,33),I=1,3),((AB(J,I,33),J=1,2),I=1,3) /
     :                             2,  0, -4,
     :                    6D-6,       -13D-6,
     :                   -8D-6,         2D-6,
     :                   14D-7,        10D-7 /
      DATA (IJSP(I,34),I=1,3),((AB(J,I,34),J=1,2),I=1,3) /
     :                             2,  0, -3,
     :                   10D-6,        22D-6,
     :                   10D-6,        -7D-6,
     :                  -65D-7,        12D-7 /
      DATA (IJSP(I,35),I=1,3),((AB(J,I,35),J=1,2),I=1,3) /
     :                             2,  0, -2,
     :                  -57D-6,       -32D-6,
     :                    0D-6,        21D-6,
     :                  126D-7,      -233D-7 /
      DATA (IJSP(I,36),I=1,3),((AB(J,I,36),J=1,2),I=1,3) /
     :                             2,  0, -1,
     :                  157D-6,       -46D-6,
     :                    8D-6,         5D-6,
     :                  270D-7,      1068D-7 /
      DATA (IJSP(I,37),I=1,3),((AB(J,I,37),J=1,2),I=1,3) /
     :                             2,  0,  0,
     :                   12D-6,       -18D-6,
     :                   13D-6,        16D-6,
     :                  254D-7,       155D-7 /
      DATA (IJSP(I,38),I=1,3),((AB(J,I,38),J=1,2),I=1,3) /
     :                             2,  0,  1,
     :                   -4D-6,         8D-6,
     :                   -2D-6,        -3D-6,
     :                  -26D-7,        -2D-7 /
      DATA (IJSP(I,39),I=1,3),((AB(J,I,39),J=1,2),I=1,3) /
     :                             2,  0,  2,
     :                   -5D-6,         0D-6,
     :                    0D-6,         0D-6,
     :                    7D-7,         0D-7 /
      DATA (IJSP(I,40),I=1,3),((AB(J,I,40),J=1,2),I=1,3) /
     :                             2,  0,  3,
     :                    3D-6,         4D-6,
     :                    0D-6,         1D-6,
     :                  -11D-7,         4D-7 /
      DATA (IJSP(I,41),I=1,3),((AB(J,I,41),J=1,2),I=1,3) /
     :                             3,  0, -2,
     :                   -1D-6,        -1D-6,
     :                    0D-6,         1D-6,
     :                    4D-7,       -14D-7 /
      DATA (IJSP(I,42),I=1,3),((AB(J,I,42),J=1,2),I=1,3) /
     :                             3,  0, -1,
     :                    6D-6,        -3D-6,
     :                    0D-6,         0D-6,
     :                   18D-7,        35D-7 /
      DATA (IJSP(I,43),I=1,3),((AB(J,I,43),J=1,2),I=1,3) /
     :                             3,  0,  0,
     :                   -1D-6,        -2D-6,
     :                    0D-6,         1D-6,
     :                   13D-7,         3D-7 /


*  Validate the planet number.
      IF (NP.LT.1.OR.NP.GT.9) THEN
         JSTAT=-1
         DO I=1,6
            PV(I)=0D0
         END DO
      ELSE

*     Separate algorithms for Pluto and the rest.
         IF (NP.NE.9) THEN

*        -----------------------
*        Mercury through Neptune
*        -----------------------

*        Time: Julian millennia since J2000.
            T=(DATE-51544.5D0)/365250D0

*        OK status unless remote epoch.
            IF (ABS(T).LE.1D0) THEN
               JSTAT=0
            ELSE
               JSTAT=1
            END IF

*        Compute the mean elements.
            DA=A(1,NP)+(A(2,NP)+A(3,NP)*T)*T
            DL=(3600D0*DLM(1,NP)+(DLM(2,NP)+DLM(3,NP)*T)*T)*AS2R
            DE=E(1,NP)+(E(2,NP)+E(3,NP)*T)*T
            DPE=MOD((3600D0*PI(1,NP)+(PI(2,NP)+PI(3,NP)*T)*T)*AS2R,D2PI)
            DI=(3600D0*DINC(1,NP)+(DINC(2,NP)+DINC(3,NP)*T)*T)*AS2R
            DO=MOD((3600D0*OMEGA(1,NP)
     :                        +(OMEGA(2,NP)+OMEGA(3,NP)*T)*T)*AS2R,D2PI)

*        Apply the trigonometric terms.
            DMU=0.35953620D0*T
            DO J=1,8
               ARGA=DKP(J,NP)*DMU
               ARGL=DKQ(J,NP)*DMU
               DA=DA+(CA(J,NP)*COS(ARGA)+SA(J,NP)*SIN(ARGA))*1D-7
               DL=DL+(CLO(J,NP)*COS(ARGL)+SLO(J,NP)*SIN(ARGL))*1D-7
            END DO
            ARGA=DKP(9,NP)*DMU
            DA=DA+T*(CA(9,NP)*COS(ARGA)+SA(9,NP)*SIN(ARGA))*1D-7
            DO J=9,10
               ARGL=DKQ(J,NP)*DMU
               DL=DL+T*(CLO(J,NP)*COS(ARGL)+SLO(J,NP)*SIN(ARGL))*1D-7
            END DO
            DL=MOD(DL,D2PI)

*        Daily motion.
            DM=GCON*SQRT((1D0+1D0/AMAS(NP))/(DA*DA*DA))

*        Make the prediction.
            CALL sla_PLANEL(DATE,1,DATE,DI,DO,DPE,DA,DE,DL,DM,PV,J)
            IF (J.LT.0) JSTAT=-2

         ELSE

*        -----
*        Pluto
*        -----

*        Time: Julian centuries since J2000.
            T=(DATE-51544.5D0)/36525D0

*        OK status unless remote epoch.
            IF (T.GE.-1.15D0.AND.T.LE.1D0) THEN
               JSTAT=0
            ELSE
               JSTAT=1
            END IF

*        Fundamental arguments (radians).
            DJ=(DJ0+DJD*T)*D2R
            DS=(DS0+DSD*T)*D2R
            DP=(DP0+DPD*T)*D2R

*        Initialize coefficients and derivatives.
            DO I=1,3
               WLBR(I)=0D0
               WLBRD(I)=0D0
            END DO

*        Term by term through Meeus Table 36.A.
            DO J=1,43

*           Argument and derivative (radians, radians per century).
               WJ=DBLE(IJSP(1,J))
               WS=DBLE(IJSP(2,J))
               WP=DBLE(IJSP(3,J))
               AL=WJ*DJ+WS*DS+WP*DP
               ALD=(WJ*DJD+WS*DSD+WP*DPD)*D2R

*           Functions of argument.
               SAL=SIN(AL)
               CAL=COS(AL)

*           Periodic terms in longitude, latitude, radius vector.
               DO I=1,3

*              A and B coefficients (deg, AU).
                  AC=AB(1,I,J)
                  BC=AB(2,I,J)

*              Periodic terms (deg, AU, deg/Jc, AU/Jc).
                  WLBR(I)=WLBR(I)+AC*SAL+BC*CAL
                  WLBRD(I)=WLBRD(I)+(AC*CAL-BC*SAL)*ALD
               END DO
            END DO

*        Heliocentric longitude and derivative (radians, radians/sec).
            DL=(DL0+DLD0*T+WLBR(1))*D2R
            DLD=(DLD0+WLBRD(1))*D2R/SPC

*        Heliocentric latitude and derivative (radians, radians/sec).
            DB=(DB0+WLBR(2))*D2R
            DBD=WLBRD(2)*D2R/SPC

*        Heliocentric radius vector and derivative (AU, AU/sec).
            DR=DR0+WLBR(3)
            DRD=WLBRD(3)/SPC

*        Functions of latitude, longitude, radius vector.
            SL=SIN(DL)
            CL=COS(DL)
            SB=SIN(DB)
            CB=COS(DB)
            SLCB=SL*CB
            CLCB=CL*CB

*        Heliocentric vector and derivative, J2000 ecliptic and equinox.
            X=DR*CLCB
            Y=DR*SLCB
            Z=DR*SB
            XD=DRD*CLCB-DR*(CL*SB*DBD+SLCB*DLD)
            YD=DRD*SLCB+DR*(-SL*SB*DBD+CLCB*DLD)
            ZD=DRD*SB+DR*CB*DBD

*        Transform to J2000 equator and equinox.
            PV(1)=X
            PV(2)=Y*CE-Z*SE
            PV(3)=Y*SE+Z*CE
            PV(4)=XD
            PV(5)=YD*CE-ZD*SE
            PV(6)=YD*SE+ZD*CE
         END IF
      END IF

      END
