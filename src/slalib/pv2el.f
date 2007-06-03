      SUBROUTINE sla_PV2EL (PV, DATE, PMASS, JFORMR,
     :                      JFORM, EPOCH, ORBINC, ANODE, PERIH,
     :                      AORQ, E, AORL, DM, JSTAT)
*+
*     - - - - - -
*      P V 2 E L
*     - - - - - -
*
*  Heliocentric osculating elements obtained from instantaneous position
*  and velocity.
*
*  Given:
*     PV        d(6)   heliocentric x,y,z,xdot,ydot,zdot of date,
*                      J2000 equatorial triad (AU,AU/s; Note 1)
*     DATE      d      date (TT Modified Julian Date = JD-2400000.5)
*     PMASS     d      mass of the planet (Sun=1; Note 2)
*     JFORMR    i      requested element set (1-3; Note 3)
*
*  Returned:
*     JFORM     d      element set actually returned (1-3; Note 4)
*     EPOCH     d      epoch of elements (TT MJD)
*     ORBINC    d      inclination (radians)
*     ANODE     d      longitude of the ascending node (radians)
*     PERIH     d      longitude or argument of perihelion (radians)
*     AORQ      d      mean distance or perihelion distance (AU)
*     E         d      eccentricity
*     AORL      d      mean anomaly or longitude (radians, JFORM=1,2 only)
*     DM        d      daily motion (radians, JFORM=1 only)
*     JSTAT     i      status:  0 = OK
*                              -1 = illegal PMASS
*                              -2 = illegal JFORMR
*                              -3 = position/velocity out of range
*
*  Notes
*
*  1  The PV 6-vector is with respect to the mean equator and equinox of
*     epoch J2000.  The orbital elements produced are with respect to
*     the J2000 ecliptic and mean equinox.
*
*  2  The mass, PMASS, is important only for the larger planets.  For
*     most purposes (e.g. asteroids) use 0D0.  Values less than zero
*     are illegal.
*
*  3  Three different element-format options are supported:
*
*     Option JFORM=1, suitable for the major planets:
*
*     EPOCH  = epoch of elements (TT MJD)
*     ORBINC = inclination i (radians)
*     ANODE  = longitude of the ascending node, big omega (radians)
*     PERIH  = longitude of perihelion, curly pi (radians)
*     AORQ   = mean distance, a (AU)
*     E      = eccentricity, e
*     AORL   = mean longitude L (radians)
*     DM     = daily motion (radians)
*
*     Option JFORM=2, suitable for minor planets:
*
*     EPOCH  = epoch of elements (TT MJD)
*     ORBINC = inclination i (radians)
*     ANODE  = longitude of the ascending node, big omega (radians)
*     PERIH  = argument of perihelion, little omega (radians)
*     AORQ   = mean distance, a (AU)
*     E      = eccentricity, e
*     AORL   = mean anomaly M (radians)
*
*     Option JFORM=3, suitable for comets:
*
*     EPOCH  = epoch of perihelion (TT MJD)
*     ORBINC = inclination i (radians)
*     ANODE  = longitude of the ascending node, big omega (radians)
*     PERIH  = argument of perihelion, little omega (radians)
*     AORQ   = perihelion distance, q (AU)
*     E      = eccentricity, e
*
*  4  It may not be possible to generate elements in the form
*     requested through JFORMR.  The caller is notified of the form
*     of elements actually returned by means of the JFORM argument:
*
*      JFORMR   JFORM     meaning
*
*        1        1       OK - elements are in the requested format
*        1        2       never happens
*        1        3       orbit not elliptical
*
*        2        1       never happens
*        2        2       OK - elements are in the requested format
*        2        3       orbit not elliptical
*
*        3        1       never happens
*        3        2       never happens
*        3        3       OK - elements are in the requested format
*
*  5  The arguments returned for each value of JFORM (cf Note 5: JFORM
*     may not be the same as JFORMR) are as follows:
*
*         JFORM         1              2              3
*         EPOCH         t0             t0             T
*         ORBINC        i              i              i
*         ANODE         Omega          Omega          Omega
*         PERIH         curly pi       omega          omega
*         AORQ          a              a              q
*         E             e              e              e
*         AORL          L              M              -
*         DM            n              -              -
*
*     where:
*
*         t0           is the epoch of the elements (MJD, TT)
*         T              "    epoch of perihelion (MJD, TT)
*         i              "    inclination (radians)
*         Omega          "    longitude of the ascending node (radians)
*         curly pi       "    longitude of perihelion (radians)
*         omega          "    argument of perihelion (radians)
*         a              "    mean distance (AU)
*         q              "    perihelion distance (AU)
*         e              "    eccentricity
*         L              "    longitude (radians, 0-2pi)
*         M              "    mean anomaly (radians, 0-2pi)
*         n              "    daily motion (radians)
*         -             means no value is set
*
*  6  At very small inclinations, the longitude of the ascending node
*     ANODE becomes indeterminate and under some circumstances may be
*     set arbitrarily to zero.  Similarly, if the orbit is close to
*     circular, the true anomaly becomes indeterminate and under some
*     circumstances may be set arbitrarily to zero.  In such cases,
*     the other elements are automatically adjusted to compensate,
*     and so the elements remain a valid description of the orbit.
*
*  7  The osculating epoch for the returned elements is the argument
*     DATE.
*
*  Reference:  Sterne, Theodore E., "An Introduction to Celestial
*              Mechanics", Interscience Publishers, 1960
*
*  Called:  sla_DRANRM
*
*  Last revision:   8 September 2005
*
*  Copyright P.T.Wallace.  All rights reserved.
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

      DOUBLE PRECISION PV(6),DATE,PMASS
      INTEGER JFORMR,JFORM
      DOUBLE PRECISION EPOCH,ORBINC,ANODE,PERIH,AORQ,E,AORL,DM
      INTEGER JSTAT

*  Seconds to days
      DOUBLE PRECISION DAY
      PARAMETER (DAY=86400D0)

*  Gaussian gravitational constant (exact)
      DOUBLE PRECISION GCON
      PARAMETER (GCON=0.01720209895D0)

*  Sin and cos of J2000 mean obliquity (IAU 1976)
      DOUBLE PRECISION SE,CE
      PARAMETER (SE=0.3977771559319137D0,
     :           CE=0.9174820620691818D0)

*  Minimum allowed distance (AU) and speed (AU/day)
      DOUBLE PRECISION RMIN,VMIN
      PARAMETER (RMIN=1D-3,VMIN=1D-8)

*  How close to unity the eccentricity has to be to call it a parabola
      DOUBLE PRECISION PARAB
      PARAMETER (PARAB=1D-8)

      DOUBLE PRECISION X,Y,Z,XD,YD,ZD,R,V2,V,RDV,GMU,HX,HY,HZ,
     :                 HX2PY2,H2,H,OI,BIGOM,AR,ECC,S,C,AT,U,OM,
     :                 GAR3,EM1,EP1,HAT,SHAT,CHAT,AE,AM,DN,PL,
     :                 EL,Q,TP,THAT,THHF,F

      INTEGER JF

      DOUBLE PRECISION sla_DRANRM


*  Validate arguments PMASS and JFORMR.
      IF (PMASS.LT.0D0) THEN
         JSTAT = -1
         GO TO 999
      END IF
      IF (JFORMR.LT.1.OR.JFORMR.GT.3) THEN
         JSTAT = -2
         GO TO 999
      END IF

*  Provisionally assume the elements will be in the chosen form.
      JF = JFORMR

*  Rotate the position from equatorial to ecliptic coordinates.
      X = PV(1)
      Y = PV(2)*CE+PV(3)*SE
      Z = -PV(2)*SE+PV(3)*CE

*  Rotate the velocity similarly, scaling to AU/day.
      XD = DAY*PV(4)
      YD = DAY*(PV(5)*CE+PV(6)*SE)
      ZD = DAY*(-PV(5)*SE+PV(6)*CE)

*  Distance and speed.
      R = SQRT(X*X+Y*Y+Z*Z)
      V2 = XD*XD+YD*YD+ZD*ZD
      V = SQRT(V2)

*  Reject unreasonably small values.
      IF (R.LT.RMIN.OR.V.LT.VMIN) THEN
         JSTAT = -3
         GO TO 999
      END IF

*  R dot V.
      RDV = X*XD+Y*YD+Z*ZD

*  Mu.
      GMU = (1D0+PMASS)*GCON*GCON

*  Vector angular momentum per unit reduced mass.
      HX = Y*ZD-Z*YD
      HY = Z*XD-X*ZD
      HZ = X*YD-Y*XD

*  Areal constant.
      HX2PY2 = HX*HX+HY*HY
      H2 = HX2PY2+HZ*HZ
      H = SQRT(H2)

*  Inclination.
      OI = ATAN2(SQRT(HX2PY2),HZ)

*  Longitude of ascending node.
      IF (HX.NE.0D0.OR.HY.NE.0D0) THEN
         BIGOM = ATAN2(HX,-HY)
      ELSE
         BIGOM=0D0
      END IF

*  Reciprocal of mean distance etc.
      AR = 2D0/R-V2/GMU

*  Eccentricity.
      ECC = SQRT(MAX(1D0-AR*H2/GMU,0D0))

*  True anomaly.
      S = H*RDV
      C = H2-R*GMU
      IF (S.NE.0D0.OR.C.NE.0D0) THEN
         AT = ATAN2(S,C)
      ELSE
         AT = 0D0
      END IF

*  Argument of the latitude.
      S = SIN(BIGOM)
      C = COS(BIGOM)
      U = ATAN2((-X*S+Y*C)*COS(OI)+Z*SIN(OI),X*C+Y*S)

*  Argument of perihelion.
      OM = U-AT

*  Capture near-parabolic cases.
      IF (ABS(ECC-1D0).LT.PARAB) ECC=1D0

*  Comply with JFORMR = 1 or 2 only if orbit is elliptical.
      IF (ECC.GE.1D0) JF=3

*  Functions.
      GAR3 = GMU*AR*AR*AR
      EM1 = ECC-1D0
      EP1 = ECC+1D0
      HAT = AT/2D0
      SHAT = SIN(HAT)
      CHAT = COS(HAT)

*  Variable initializations to avoid compiler warnings.
      AM = 0D0
      DN = 0D0
      PL = 0D0
      EL = 0D0
      Q = 0D0
      TP = 0D0

*  Ellipse?
      IF (ECC.LT.1D0 ) THEN

*     Eccentric anomaly.
         AE = 2D0*ATAN2(SQRT(-EM1)*SHAT,SQRT(EP1)*CHAT)

*     Mean anomaly.
         AM = AE-ECC*SIN(AE)

*     Daily motion.
         DN = SQRT(GAR3)
      END IF

*  "Major planet" element set?
      IF (JF.EQ.1) THEN

*     Longitude of perihelion.
         PL = BIGOM+OM

*     Longitude at epoch.
         EL = PL+AM
      END IF

*  "Comet" element set?
      IF (JF.EQ.3) THEN

*     Perihelion distance.
         Q = H2/(GMU*EP1)

*     Ellipse, parabola, hyperbola?
         IF (ECC.LT.1D0) THEN

*        Ellipse: epoch of perihelion.
            TP = DATE-AM/DN
         ELSE

*        Parabola or hyperbola: evaluate tan ( ( true anomaly ) / 2 )
            THAT = SHAT/CHAT
            IF (ECC.EQ.1D0) THEN

*           Parabola: epoch of perihelion.
               TP = DATE-THAT*(1D0+THAT*THAT/3D0)*H*H2/(2D0*GMU*GMU)
            ELSE

*           Hyperbola: epoch of perihelion.
               THHF = SQRT(EM1/EP1)*THAT
               F = LOG(1D0+THHF)-LOG(1D0-THHF)
               TP = DATE-(ECC*SINH(F)-F)/SQRT(-GAR3)
            END IF
         END IF
      END IF

*  Return the appropriate set of elements.
      JFORM = JF
      ORBINC = OI
      ANODE = sla_DRANRM(BIGOM)
      E = ECC
      IF (JF.EQ.1) THEN
         PERIH = sla_DRANRM(PL)
         AORL = sla_DRANRM(EL)
         DM = DN
      ELSE
         PERIH = sla_DRANRM(OM)
         IF (JF.EQ.2) AORL = sla_DRANRM(AM)
      END IF
      IF (JF.NE.3) THEN
         EPOCH = DATE
         AORQ = 1D0/AR
      ELSE
         EPOCH = TP
         AORQ = Q
      END IF
      JSTAT = 0

 999  CONTINUE
      END
