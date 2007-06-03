      SUBROUTINE sla_DMOON (DATE, PV)
*+
*     - - - - - -
*      D M O O N
*     - - - - - -
*
*  Approximate geocentric position and velocity of the Moon
*  (double precision)
*
*  Given:
*     DATE       D       TDB (loosely ET) as a Modified Julian Date
*                                                    (JD-2400000.5)
*
*  Returned:
*     PV         D(6)    Moon x,y,z,xdot,ydot,zdot, mean equator and
*                                         equinox of date (AU, AU/s)
*
*  Notes:
*
*  1  This routine is a full implementation of the algorithm
*     published by Meeus (see reference).
*
*  2  Meeus quotes accuracies of 10 arcsec in longitude, 3 arcsec in
*     latitude and 0.2 arcsec in HP (equivalent to about 20 km in
*     distance).  Comparison with JPL DE200 over the interval
*     1960-2025 gives RMS errors of 3.7 arcsec and 83 mas/hour in
*     longitude, 2.3 arcsec and 48 mas/hour in latitude, 11 km
*     and 81 mm/s in distance.  The maximum errors over the same
*     interval are 18 arcsec and 0.50 arcsec/hour in longitude,
*     11 arcsec and 0.24 arcsec/hour in latitude, 40 km and 0.29 m/s
*     in distance.
*
*  3  The original algorithm is expressed in terms of the obsolete
*     timescale Ephemeris Time.  Either TDB or TT can be used, but
*     not UT without incurring significant errors (30 arcsec at
*     the present time) due to the Moon's 0.5 arcsec/sec movement.
*
*  4  The algorithm is based on pre IAU 1976 standards.  However,
*     the result has been moved onto the new (FK5) equinox, an
*     adjustment which is in any case much smaller than the
*     intrinsic accuracy of the procedure.
*
*  5  Velocity is obtained by a complete analytical differentiation
*     of the Meeus model.
*
*  Reference:
*     Meeus, l'Astronomie, June 1984, p348.
*
*  P.T.Wallace   Starlink   22 January 1998
*
*  Copyright (C) 1998 Rutherford Appleton Laboratory
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

      DOUBLE PRECISION DATE,PV(6)

*  Degrees, arcseconds and seconds of time to radians
      DOUBLE PRECISION D2R,DAS2R,DS2R
      PARAMETER (D2R=0.0174532925199432957692369D0,
     :           DAS2R=4.848136811095359935899141D-6,
     :           DS2R=7.272205216643039903848712D-5)

*  Seconds per Julian century (86400*36525)
      DOUBLE PRECISION CJ
      PARAMETER (CJ=3155760000D0)

*  Julian epoch of B1950
      DOUBLE PRECISION B1950
      PARAMETER (B1950=1949.9997904423D0)

*  Earth equatorial radius in AU ( = 6378.137 / 149597870 )
      DOUBLE PRECISION ERADAU
      PARAMETER (ERADAU=4.2635212653763D-5)

      DOUBLE PRECISION T,THETA,SINOM,COSOM,DOMCOM,WA,DWA,WB,DWB,WOM,
     :                 DWOM,SINWOM,COSWOM,V,DV,COEFF,EMN,EMPN,DN,FN,EN,
     :                 DEN,DTHETA,FTHETA,EL,DEL,B,DB,BF,DBF,P,DP,SP,R,
     :                 DR,X,Y,Z,XD,YD,ZD,SEL,CEL,SB,CB,RCB,RBD,W,EPJ,
     :                 EQCOR,EPS,SINEPS,COSEPS,ES,EC
      INTEGER N,I

*
*  Coefficients for fundamental arguments
*
*   at J1900:  T**0, T**1, T**2, T**3
*   at epoch:  T**0, T**1
*
*  Units are degrees for position and Julian centuries for time
*

*  Moon's mean longitude
      DOUBLE PRECISION ELP0,ELP1,ELP2,ELP3,ELP,DELP
      PARAMETER (ELP0=270.434164D0,
     :           ELP1=481267.8831D0,
     :           ELP2=-0.001133D0,
     :           ELP3=0.0000019D0)

*  Sun's mean anomaly
      DOUBLE PRECISION EM0,EM1,EM2,EM3,EM,DEM
      PARAMETER (EM0=358.475833D0,
     :           EM1=35999.0498D0,
     :           EM2=-0.000150D0,
     :           EM3=-0.0000033D0)

*  Moon's mean anomaly
      DOUBLE PRECISION EMP0,EMP1,EMP2,EMP3,EMP,DEMP
      PARAMETER (EMP0=296.104608D0,
     :           EMP1=477198.8491D0,
     :           EMP2=0.009192D0,
     :           EMP3=0.0000144D0)

*  Moon's mean elongation
      DOUBLE PRECISION D0,D1,D2,D3,D,DD
      PARAMETER (D0=350.737486D0,
     :           D1=445267.1142D0,
     :           D2=-0.001436D0,
     :           D3=0.0000019D0)

*  Mean distance of the Moon from its ascending node
      DOUBLE PRECISION F0,F1,F2,F3,F,DF
      PARAMETER (F0=11.250889D0,
     :           F1=483202.0251D0,
     :           F2=-0.003211D0,
     :           F3=-0.0000003D0)

*  Longitude of the Moon's ascending node
      DOUBLE PRECISION OM0,OM1,OM2,OM3,OM,DOM
      PARAMETER (OM0=259.183275D0,
     :           OM1=-1934.1420D0,
     :           OM2=0.002078D0,
     :           OM3=0.0000022D0)

*  Coefficients for (dimensionless) E factor
      DOUBLE PRECISION E1,E2,E,DE,ESQ,DESQ
      PARAMETER (E1=-0.002495D0,E2=-0.00000752D0)

*  Coefficients for periodic variations etc
      DOUBLE PRECISION PAC,PA0,PA1
      PARAMETER (PAC=0.000233D0,PA0=51.2D0,PA1=20.2D0)
      DOUBLE PRECISION PBC
      PARAMETER (PBC=-0.001778D0)
      DOUBLE PRECISION PCC
      PARAMETER (PCC=0.000817D0)
      DOUBLE PRECISION PDC
      PARAMETER (PDC=0.002011D0)
      DOUBLE PRECISION PEC,PE0,PE1,PE2
      PARAMETER (PEC=0.003964D0,
     :                     PE0=346.560D0,PE1=132.870D0,PE2=-0.0091731D0)
      DOUBLE PRECISION PFC
      PARAMETER (PFC=0.001964D0)
      DOUBLE PRECISION PGC
      PARAMETER (PGC=0.002541D0)
      DOUBLE PRECISION PHC
      PARAMETER (PHC=0.001964D0)
      DOUBLE PRECISION PIC
      PARAMETER (PIC=-0.024691D0)
      DOUBLE PRECISION PJC,PJ0,PJ1
      PARAMETER (PJC=-0.004328D0,PJ0=275.05D0,PJ1=-2.30D0)
      DOUBLE PRECISION CW1
      PARAMETER (CW1=0.0004664D0)
      DOUBLE PRECISION CW2
      PARAMETER (CW2=0.0000754D0)

*
*  Coefficients for Moon position
*
*   Tx(N)       = coefficient of L, B or P term (deg)
*   ITx(N,1-5)  = coefficients of M, M', D, F, E**n in argument
*
      INTEGER NL,NB,NP
      PARAMETER (NL=50,NB=45,NP=31)
      DOUBLE PRECISION TL(NL),TB(NB),TP(NP)
      INTEGER ITL(5,NL),ITB(5,NB),ITP(5,NP)
*
*  Longitude
*                                         M   M'  D   F   n
      DATA TL( 1)/            +6.288750D0                     /,
     :     (ITL(I, 1),I=1,5)/            +0, +1, +0, +0,  0   /
      DATA TL( 2)/            +1.274018D0                     /,
     :     (ITL(I, 2),I=1,5)/            +0, -1, +2, +0,  0   /
      DATA TL( 3)/            +0.658309D0                     /,
     :     (ITL(I, 3),I=1,5)/            +0, +0, +2, +0,  0   /
      DATA TL( 4)/            +0.213616D0                     /,
     :     (ITL(I, 4),I=1,5)/            +0, +2, +0, +0,  0   /
      DATA TL( 5)/            -0.185596D0                     /,
     :     (ITL(I, 5),I=1,5)/            +1, +0, +0, +0,  1   /
      DATA TL( 6)/            -0.114336D0                     /,
     :     (ITL(I, 6),I=1,5)/            +0, +0, +0, +2,  0   /
      DATA TL( 7)/            +0.058793D0                     /,
     :     (ITL(I, 7),I=1,5)/            +0, -2, +2, +0,  0   /
      DATA TL( 8)/            +0.057212D0                     /,
     :     (ITL(I, 8),I=1,5)/            -1, -1, +2, +0,  1   /
      DATA TL( 9)/            +0.053320D0                     /,
     :     (ITL(I, 9),I=1,5)/            +0, +1, +2, +0,  0   /
      DATA TL(10)/            +0.045874D0                     /,
     :     (ITL(I,10),I=1,5)/            -1, +0, +2, +0,  1   /
      DATA TL(11)/            +0.041024D0                     /,
     :     (ITL(I,11),I=1,5)/            -1, +1, +0, +0,  1   /
      DATA TL(12)/            -0.034718D0                     /,
     :     (ITL(I,12),I=1,5)/            +0, +0, +1, +0,  0   /
      DATA TL(13)/            -0.030465D0                     /,
     :     (ITL(I,13),I=1,5)/            +1, +1, +0, +0,  1   /
      DATA TL(14)/            +0.015326D0                     /,
     :     (ITL(I,14),I=1,5)/            +0, +0, +2, -2,  0   /
      DATA TL(15)/            -0.012528D0                     /,
     :     (ITL(I,15),I=1,5)/            +0, +1, +0, +2,  0   /
      DATA TL(16)/            -0.010980D0                     /,
     :     (ITL(I,16),I=1,5)/            +0, -1, +0, +2,  0   /
      DATA TL(17)/            +0.010674D0                     /,
     :     (ITL(I,17),I=1,5)/            +0, -1, +4, +0,  0   /
      DATA TL(18)/            +0.010034D0                     /,
     :     (ITL(I,18),I=1,5)/            +0, +3, +0, +0,  0   /
      DATA TL(19)/            +0.008548D0                     /,
     :     (ITL(I,19),I=1,5)/            +0, -2, +4, +0,  0   /
      DATA TL(20)/            -0.007910D0                     /,
     :     (ITL(I,20),I=1,5)/            +1, -1, +2, +0,  1   /
      DATA TL(21)/            -0.006783D0                     /,
     :     (ITL(I,21),I=1,5)/            +1, +0, +2, +0,  1   /
      DATA TL(22)/            +0.005162D0                     /,
     :     (ITL(I,22),I=1,5)/            +0, +1, -1, +0,  0   /
      DATA TL(23)/            +0.005000D0                     /,
     :     (ITL(I,23),I=1,5)/            +1, +0, +1, +0,  1   /
      DATA TL(24)/            +0.004049D0                     /,
     :     (ITL(I,24),I=1,5)/            -1, +1, +2, +0,  1   /
      DATA TL(25)/            +0.003996D0                     /,
     :     (ITL(I,25),I=1,5)/            +0, +2, +2, +0,  0   /
      DATA TL(26)/            +0.003862D0                     /,
     :     (ITL(I,26),I=1,5)/            +0, +0, +4, +0,  0   /
      DATA TL(27)/            +0.003665D0                     /,
     :     (ITL(I,27),I=1,5)/            +0, -3, +2, +0,  0   /
      DATA TL(28)/            +0.002695D0                     /,
     :     (ITL(I,28),I=1,5)/            -1, +2, +0, +0,  1   /
      DATA TL(29)/            +0.002602D0                     /,
     :     (ITL(I,29),I=1,5)/            +0, +1, -2, -2,  0   /
      DATA TL(30)/            +0.002396D0                     /,
     :     (ITL(I,30),I=1,5)/            -1, -2, +2, +0,  1   /
      DATA TL(31)/            -0.002349D0                     /,
     :     (ITL(I,31),I=1,5)/            +0, +1, +1, +0,  0   /
      DATA TL(32)/            +0.002249D0                     /,
     :     (ITL(I,32),I=1,5)/            -2, +0, +2, +0,  2   /
      DATA TL(33)/            -0.002125D0                     /,
     :     (ITL(I,33),I=1,5)/            +1, +2, +0, +0,  1   /
      DATA TL(34)/            -0.002079D0                     /,
     :     (ITL(I,34),I=1,5)/            +2, +0, +0, +0,  2   /
      DATA TL(35)/            +0.002059D0                     /,
     :     (ITL(I,35),I=1,5)/            -2, -1, +2, +0,  2   /
      DATA TL(36)/            -0.001773D0                     /,
     :     (ITL(I,36),I=1,5)/            +0, +1, +2, -2,  0   /
      DATA TL(37)/            -0.001595D0                     /,
     :     (ITL(I,37),I=1,5)/            +0, +0, +2, +2,  0   /
      DATA TL(38)/            +0.001220D0                     /,
     :     (ITL(I,38),I=1,5)/            -1, -1, +4, +0,  1   /
      DATA TL(39)/            -0.001110D0                     /,
     :     (ITL(I,39),I=1,5)/            +0, +2, +0, +2,  0   /
      DATA TL(40)/            +0.000892D0                     /,
     :     (ITL(I,40),I=1,5)/            +0, +1, -3, +0,  0   /
      DATA TL(41)/            -0.000811D0                     /,
     :     (ITL(I,41),I=1,5)/            +1, +1, +2, +0,  1   /
      DATA TL(42)/            +0.000761D0                     /,
     :     (ITL(I,42),I=1,5)/            -1, -2, +4, +0,  1   /
      DATA TL(43)/            +0.000717D0                     /,
     :     (ITL(I,43),I=1,5)/            -2, +1, +0, +0,  2   /
      DATA TL(44)/            +0.000704D0                     /,
     :     (ITL(I,44),I=1,5)/            -2, +1, -2, +0,  2   /
      DATA TL(45)/            +0.000693D0                     /,
     :     (ITL(I,45),I=1,5)/            +1, -2, +2, +0,  1   /
      DATA TL(46)/            +0.000598D0                     /,
     :     (ITL(I,46),I=1,5)/            -1, +0, +2, -2,  1   /
      DATA TL(47)/            +0.000550D0                     /,
     :     (ITL(I,47),I=1,5)/            +0, +1, +4, +0,  0   /
      DATA TL(48)/            +0.000538D0                     /,
     :     (ITL(I,48),I=1,5)/            +0, +4, +0, +0,  0   /
      DATA TL(49)/            +0.000521D0                     /,
     :     (ITL(I,49),I=1,5)/            -1, +0, +4, +0,  1   /
      DATA TL(50)/            +0.000486D0                     /,
     :     (ITL(I,50),I=1,5)/            +0, +2, -1, +0,  0   /
*
*  Latitude
*                                         M   M'  D   F   n
      DATA TB( 1)/            +5.128189D0                     /,
     :     (ITB(I, 1),I=1,5)/            +0, +0, +0, +1,  0   /
      DATA TB( 2)/            +0.280606D0                     /,
     :     (ITB(I, 2),I=1,5)/            +0, +1, +0, +1,  0   /
      DATA TB( 3)/            +0.277693D0                     /,
     :     (ITB(I, 3),I=1,5)/            +0, +1, +0, -1,  0   /
      DATA TB( 4)/            +0.173238D0                     /,
     :     (ITB(I, 4),I=1,5)/            +0, +0, +2, -1,  0   /
      DATA TB( 5)/            +0.055413D0                     /,
     :     (ITB(I, 5),I=1,5)/            +0, -1, +2, +1,  0   /
      DATA TB( 6)/            +0.046272D0                     /,
     :     (ITB(I, 6),I=1,5)/            +0, -1, +2, -1,  0   /
      DATA TB( 7)/            +0.032573D0                     /,
     :     (ITB(I, 7),I=1,5)/            +0, +0, +2, +1,  0   /
      DATA TB( 8)/            +0.017198D0                     /,
     :     (ITB(I, 8),I=1,5)/            +0, +2, +0, +1,  0   /
      DATA TB( 9)/            +0.009267D0                     /,
     :     (ITB(I, 9),I=1,5)/            +0, +1, +2, -1,  0   /
      DATA TB(10)/            +0.008823D0                     /,
     :     (ITB(I,10),I=1,5)/            +0, +2, +0, -1,  0   /
      DATA TB(11)/            +0.008247D0                     /,
     :     (ITB(I,11),I=1,5)/            -1, +0, +2, -1,  1   /
      DATA TB(12)/            +0.004323D0                     /,
     :     (ITB(I,12),I=1,5)/            +0, -2, +2, -1,  0   /
      DATA TB(13)/            +0.004200D0                     /,
     :     (ITB(I,13),I=1,5)/            +0, +1, +2, +1,  0   /
      DATA TB(14)/            +0.003372D0                     /,
     :     (ITB(I,14),I=1,5)/            -1, +0, -2, +1,  1   /
      DATA TB(15)/            +0.002472D0                     /,
     :     (ITB(I,15),I=1,5)/            -1, -1, +2, +1,  1   /
      DATA TB(16)/            +0.002222D0                     /,
     :     (ITB(I,16),I=1,5)/            -1, +0, +2, +1,  1   /
      DATA TB(17)/            +0.002072D0                     /,
     :     (ITB(I,17),I=1,5)/            -1, -1, +2, -1,  1   /
      DATA TB(18)/            +0.001877D0                     /,
     :     (ITB(I,18),I=1,5)/            -1, +1, +0, +1,  1   /
      DATA TB(19)/            +0.001828D0                     /,
     :     (ITB(I,19),I=1,5)/            +0, -1, +4, -1,  0   /
      DATA TB(20)/            -0.001803D0                     /,
     :     (ITB(I,20),I=1,5)/            +1, +0, +0, +1,  1   /
      DATA TB(21)/            -0.001750D0                     /,
     :     (ITB(I,21),I=1,5)/            +0, +0, +0, +3,  0   /
      DATA TB(22)/            +0.001570D0                     /,
     :     (ITB(I,22),I=1,5)/            -1, +1, +0, -1,  1   /
      DATA TB(23)/            -0.001487D0                     /,
     :     (ITB(I,23),I=1,5)/            +0, +0, +1, +1,  0   /
      DATA TB(24)/            -0.001481D0                     /,
     :     (ITB(I,24),I=1,5)/            +1, +1, +0, +1,  1   /
      DATA TB(25)/            +0.001417D0                     /,
     :     (ITB(I,25),I=1,5)/            -1, -1, +0, +1,  1   /
      DATA TB(26)/            +0.001350D0                     /,
     :     (ITB(I,26),I=1,5)/            -1, +0, +0, +1,  1   /
      DATA TB(27)/            +0.001330D0                     /,
     :     (ITB(I,27),I=1,5)/            +0, +0, -1, +1,  0   /
      DATA TB(28)/            +0.001106D0                     /,
     :     (ITB(I,28),I=1,5)/            +0, +3, +0, +1,  0   /
      DATA TB(29)/            +0.001020D0                     /,
     :     (ITB(I,29),I=1,5)/            +0, +0, +4, -1,  0   /
      DATA TB(30)/            +0.000833D0                     /,
     :     (ITB(I,30),I=1,5)/            +0, -1, +4, +1,  0   /
      DATA TB(31)/            +0.000781D0                     /,
     :     (ITB(I,31),I=1,5)/            +0, +1, +0, -3,  0   /
      DATA TB(32)/            +0.000670D0                     /,
     :     (ITB(I,32),I=1,5)/            +0, -2, +4, +1,  0   /
      DATA TB(33)/            +0.000606D0                     /,
     :     (ITB(I,33),I=1,5)/            +0, +0, +2, -3,  0   /
      DATA TB(34)/            +0.000597D0                     /,
     :     (ITB(I,34),I=1,5)/            +0, +2, +2, -1,  0   /
      DATA TB(35)/            +0.000492D0                     /,
     :     (ITB(I,35),I=1,5)/            -1, +1, +2, -1,  1   /
      DATA TB(36)/            +0.000450D0                     /,
     :     (ITB(I,36),I=1,5)/            +0, +2, -2, -1,  0   /
      DATA TB(37)/            +0.000439D0                     /,
     :     (ITB(I,37),I=1,5)/            +0, +3, +0, -1,  0   /
      DATA TB(38)/            +0.000423D0                     /,
     :     (ITB(I,38),I=1,5)/            +0, +2, +2, +1,  0   /
      DATA TB(39)/            +0.000422D0                     /,
     :     (ITB(I,39),I=1,5)/            +0, -3, +2, -1,  0   /
      DATA TB(40)/            -0.000367D0                     /,
     :     (ITB(I,40),I=1,5)/            +1, -1, +2, +1,  1   /
      DATA TB(41)/            -0.000353D0                     /,
     :     (ITB(I,41),I=1,5)/            +1, +0, +2, +1,  1   /
      DATA TB(42)/            +0.000331D0                     /,
     :     (ITB(I,42),I=1,5)/            +0, +0, +4, +1,  0   /
      DATA TB(43)/            +0.000317D0                     /,
     :     (ITB(I,43),I=1,5)/            -1, +1, +2, +1,  1   /
      DATA TB(44)/            +0.000306D0                     /,
     :     (ITB(I,44),I=1,5)/            -2, +0, +2, -1,  2   /
      DATA TB(45)/            -0.000283D0                     /,
     :     (ITB(I,45),I=1,5)/            +0, +1, +0, +3,  0   /
*
*  Parallax
*                                         M   M'  D   F   n
      DATA TP( 1)/            +0.950724D0                     /,
     :     (ITP(I, 1),I=1,5)/            +0, +0, +0, +0,  0   /
      DATA TP( 2)/            +0.051818D0                     /,
     :     (ITP(I, 2),I=1,5)/            +0, +1, +0, +0,  0   /
      DATA TP( 3)/            +0.009531D0                     /,
     :     (ITP(I, 3),I=1,5)/            +0, -1, +2, +0,  0   /
      DATA TP( 4)/            +0.007843D0                     /,
     :     (ITP(I, 4),I=1,5)/            +0, +0, +2, +0,  0   /
      DATA TP( 5)/            +0.002824D0                     /,
     :     (ITP(I, 5),I=1,5)/            +0, +2, +0, +0,  0   /
      DATA TP( 6)/            +0.000857D0                     /,
     :     (ITP(I, 6),I=1,5)/            +0, +1, +2, +0,  0   /
      DATA TP( 7)/            +0.000533D0                     /,
     :     (ITP(I, 7),I=1,5)/            -1, +0, +2, +0,  1   /
      DATA TP( 8)/            +0.000401D0                     /,
     :     (ITP(I, 8),I=1,5)/            -1, -1, +2, +0,  1   /
      DATA TP( 9)/            +0.000320D0                     /,
     :     (ITP(I, 9),I=1,5)/            -1, +1, +0, +0,  1   /
      DATA TP(10)/            -0.000271D0                     /,
     :     (ITP(I,10),I=1,5)/            +0, +0, +1, +0,  0   /
      DATA TP(11)/            -0.000264D0                     /,
     :     (ITP(I,11),I=1,5)/            +1, +1, +0, +0,  1   /
      DATA TP(12)/            -0.000198D0                     /,
     :     (ITP(I,12),I=1,5)/            +0, -1, +0, +2,  0   /
      DATA TP(13)/            +0.000173D0                     /,
     :     (ITP(I,13),I=1,5)/            +0, +3, +0, +0,  0   /
      DATA TP(14)/            +0.000167D0                     /,
     :     (ITP(I,14),I=1,5)/            +0, -1, +4, +0,  0   /
      DATA TP(15)/            -0.000111D0                     /,
     :     (ITP(I,15),I=1,5)/            +1, +0, +0, +0,  1   /
      DATA TP(16)/            +0.000103D0                     /,
     :     (ITP(I,16),I=1,5)/            +0, -2, +4, +0,  0   /
      DATA TP(17)/            -0.000084D0                     /,
     :     (ITP(I,17),I=1,5)/            +0, +2, -2, +0,  0   /
      DATA TP(18)/            -0.000083D0                     /,
     :     (ITP(I,18),I=1,5)/            +1, +0, +2, +0,  1   /
      DATA TP(19)/            +0.000079D0                     /,
     :     (ITP(I,19),I=1,5)/            +0, +2, +2, +0,  0   /
      DATA TP(20)/            +0.000072D0                     /,
     :     (ITP(I,20),I=1,5)/            +0, +0, +4, +0,  0   /
      DATA TP(21)/            +0.000064D0                     /,
     :     (ITP(I,21),I=1,5)/            -1, +1, +2, +0,  1   /
      DATA TP(22)/            -0.000063D0                     /,
     :     (ITP(I,22),I=1,5)/            +1, -1, +2, +0,  1   /
      DATA TP(23)/            +0.000041D0                     /,
     :     (ITP(I,23),I=1,5)/            +1, +0, +1, +0,  1   /
      DATA TP(24)/            +0.000035D0                     /,
     :     (ITP(I,24),I=1,5)/            -1, +2, +0, +0,  1   /
      DATA TP(25)/            -0.000033D0                     /,
     :     (ITP(I,25),I=1,5)/            +0, +3, -2, +0,  0   /
      DATA TP(26)/            -0.000030D0                     /,
     :     (ITP(I,26),I=1,5)/            +0, +1, +1, +0,  0   /
      DATA TP(27)/            -0.000029D0                     /,
     :     (ITP(I,27),I=1,5)/            +0, +0, -2, +2,  0   /
      DATA TP(28)/            -0.000029D0                     /,
     :     (ITP(I,28),I=1,5)/            +1, +2, +0, +0,  1   /
      DATA TP(29)/            +0.000026D0                     /,
     :     (ITP(I,29),I=1,5)/            -2, +0, +2, +0,  2   /
      DATA TP(30)/            -0.000023D0                     /,
     :     (ITP(I,30),I=1,5)/            +0, +1, -2, +2,  0   /
      DATA TP(31)/            +0.000019D0                     /,
     :     (ITP(I,31),I=1,5)/            -1, -1, +4, +0,  1   /



*  Centuries since J1900
      T=(DATE-15019.5D0)/36525D0

*
*  Fundamental arguments (radians) and derivatives (radians per
*  Julian century) for the current epoch
*

*  Moon's mean longitude
      ELP=D2R*MOD(ELP0+(ELP1+(ELP2+ELP3*T)*T)*T,360D0)
      DELP=D2R*(ELP1+(2D0*ELP2+3D0*ELP3*T)*T)

*  Sun's mean anomaly
      EM=D2R*MOD(EM0+(EM1+(EM2+EM3*T)*T)*T,360D0)
      DEM=D2R*(EM1+(2D0*EM2+3D0*EM3*T)*T)

*  Moon's mean anomaly
      EMP=D2R*MOD(EMP0+(EMP1+(EMP2+EMP3*T)*T)*T,360D0)
      DEMP=D2R*(EMP1+(2D0*EMP2+3D0*EMP3*T)*T)

*  Moon's mean elongation
      D=D2R*MOD(D0+(D1+(D2+D3*T)*T)*T,360D0)
      DD=D2R*(D1+(2D0*D2+3D0*D3*T)*T)

*  Mean distance of the Moon from its ascending node
      F=D2R*MOD(F0+(F1+(F2+F3*T)*T)*T,360D0)
      DF=D2R*(F1+(2D0*F2+3D0*F3*T)*T)

*  Longitude of the Moon's ascending node
      OM=D2R*MOD(OM0+(OM1+(OM2+OM3*T)*T)*T,360D0)
      DOM=D2R*(OM1+(2D0*OM2+3D0*OM3*T)*T)
      SINOM=SIN(OM)
      COSOM=COS(OM)
      DOMCOM=DOM*COSOM

*  Add the periodic variations
      THETA=D2R*(PA0+PA1*T)
      WA=SIN(THETA)
      DWA=D2R*PA1*COS(THETA)
      THETA=D2R*(PE0+(PE1+PE2*T)*T)
      WB=PEC*SIN(THETA)
      DWB=D2R*PEC*(PE1+2D0*PE2*T)*COS(THETA)
      ELP=ELP+D2R*(PAC*WA+WB+PFC*SINOM)
      DELP=DELP+D2R*(PAC*DWA+DWB+PFC*DOMCOM)
      EM=EM+D2R*PBC*WA
      DEM=DEM+D2R*PBC*DWA
      EMP=EMP+D2R*(PCC*WA+WB+PGC*SINOM)
      DEMP=DEMP+D2R*(PCC*DWA+DWB+PGC*DOMCOM)
      D=D+D2R*(PDC*WA+WB+PHC*SINOM)
      DD=DD+D2R*(PDC*DWA+DWB+PHC*DOMCOM)
      WOM=OM+D2R*(PJ0+PJ1*T)
      DWOM=DOM+D2R*PJ1
      SINWOM=SIN(WOM)
      COSWOM=COS(WOM)
      F=F+D2R*(WB+PIC*SINOM+PJC*SINWOM)
      DF=DF+D2R*(DWB+PIC*DOMCOM+PJC*DWOM*COSWOM)

*  E-factor, and square
      E=1D0+(E1+E2*T)*T
      DE=E1+2D0*E2*T
      ESQ=E*E
      DESQ=2D0*E*DE

*
*  Series expansions
*

*  Longitude
      V=0D0
      DV=0D0
      DO N=NL,1,-1
         COEFF=TL(N)
         EMN=DBLE(ITL(1,N))
         EMPN=DBLE(ITL(2,N))
         DN=DBLE(ITL(3,N))
         FN=DBLE(ITL(4,N))
         I=ITL(5,N)
         IF (I.EQ.0) THEN
            EN=1D0
            DEN=0D0
         ELSE IF (I.EQ.1) THEN
            EN=E
            DEN=DE
         ELSE
            EN=ESQ
            DEN=DESQ
         END IF
         THETA=EMN*EM+EMPN*EMP+DN*D+FN*F
         DTHETA=EMN*DEM+EMPN*DEMP+DN*DD+FN*DF
         FTHETA=SIN(THETA)
         V=V+COEFF*FTHETA*EN
         DV=DV+COEFF*(COS(THETA)*DTHETA*EN+FTHETA*DEN)
      END DO
      EL=ELP+D2R*V
      DEL=(DELP+D2R*DV)/CJ

*  Latitude
      V=0D0
      DV=0D0
      DO N=NB,1,-1
         COEFF=TB(N)
         EMN=DBLE(ITB(1,N))
         EMPN=DBLE(ITB(2,N))
         DN=DBLE(ITB(3,N))
         FN=DBLE(ITB(4,N))
         I=ITB(5,N)
         IF (I.EQ.0) THEN
            EN=1D0
            DEN=0D0
         ELSE IF (I.EQ.1) THEN
            EN=E
            DEN=DE
         ELSE
            EN=ESQ
            DEN=DESQ
         END IF
         THETA=EMN*EM+EMPN*EMP+DN*D+FN*F
         DTHETA=EMN*DEM+EMPN*DEMP+DN*DD+FN*DF
         FTHETA=SIN(THETA)
         V=V+COEFF*FTHETA*EN
         DV=DV+COEFF*(COS(THETA)*DTHETA*EN+FTHETA*DEN)
      END DO
      BF=1D0-CW1*COSOM-CW2*COSWOM
      DBF=CW1*DOM*SINOM+CW2*DWOM*SINWOM
      B=D2R*V*BF
      DB=D2R*(DV*BF+V*DBF)/CJ

*  Parallax
      V=0D0
      DV=0D0
      DO N=NP,1,-1
         COEFF=TP(N)
         EMN=DBLE(ITP(1,N))
         EMPN=DBLE(ITP(2,N))
         DN=DBLE(ITP(3,N))
         FN=DBLE(ITP(4,N))
         I=ITP(5,N)
         IF (I.EQ.0) THEN
            EN=1D0
            DEN=0D0
         ELSE IF (I.EQ.1) THEN
            EN=E
            DEN=DE
         ELSE
            EN=ESQ
            DEN=DESQ
         END IF
         THETA=EMN*EM+EMPN*EMP+DN*D+FN*F
         DTHETA=EMN*DEM+EMPN*DEMP+DN*DD+FN*DF
         FTHETA=COS(THETA)
         V=V+COEFF*FTHETA*EN
         DV=DV+COEFF*(-SIN(THETA)*DTHETA*EN+FTHETA*DEN)
      END DO
      P=D2R*V
      DP=D2R*DV/CJ

*
*  Transformation into final form
*

*  Parallax to distance (AU, AU/sec)
      SP=SIN(P)
      R=ERADAU/SP
      DR=-R*DP*COS(P)/SP

*  Longitude, latitude to x,y,z (AU)
      SEL=SIN(EL)
      CEL=COS(EL)
      SB=SIN(B)
      CB=COS(B)
      RCB=R*CB
      RBD=R*DB
      W=RBD*SB-CB*DR
      X=RCB*CEL
      Y=RCB*SEL
      Z=R*SB
      XD=-Y*DEL-W*CEL
      YD=X*DEL-W*SEL
      ZD=RBD*CB+SB*DR

*  Julian centuries since J2000
      T=(DATE-51544.5D0)/36525D0

*  Fricke equinox correction
      EPJ=2000D0+T*100D0
      EQCOR=DS2R*(0.035D0+0.00085D0*(EPJ-B1950))

*  Mean obliquity (IAU 1976)
      EPS=DAS2R*(84381.448D0+(-46.8150D0+(-0.00059D0+0.001813D0*T)*T)*T)

*  To the equatorial system, mean of date, FK5 system
      SINEPS=SIN(EPS)
      COSEPS=COS(EPS)
      ES=EQCOR*SINEPS
      EC=EQCOR*COSEPS
      PV(1)=X-EC*Y+ES*Z
      PV(2)=EQCOR*X+Y*COSEPS-Z*SINEPS
      PV(3)=Y*SINEPS+Z*COSEPS
      PV(4)=XD-EC*YD+ES*ZD
      PV(5)=EQCOR*XD+YD*COSEPS-ZD*SINEPS
      PV(6)=YD*SINEPS+ZD*COSEPS

      END
