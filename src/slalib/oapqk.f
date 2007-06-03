      SUBROUTINE sla_OAPQK (TYPE, OB1, OB2, AOPRMS, RAP, DAP)
*+
*     - - - - - -
*      O A P Q K
*     - - - - - -
*
*  Quick observed to apparent place
*
*  Given:
*     TYPE   c*(*)  type of coordinates - 'R', 'H' or 'A' (see below)
*     OB1    d      observed Az, HA or RA (radians; Az is N=0,E=90)
*     OB2    d      observed ZD or Dec (radians)
*     AOPRMS d(14)  star-independent apparent-to-observed parameters:
*
*       (1)      geodetic latitude (radians)
*       (2,3)    sine and cosine of geodetic latitude
*       (4)      magnitude of diurnal aberration vector
*       (5)      height (HM)
*       (6)      ambient temperature (T)
*       (7)      pressure (P)
*       (8)      relative humidity (RH)
*       (9)      wavelength (WL)
*       (10)     lapse rate (TLR)
*       (11,12)  refraction constants A and B (radians)
*       (13)     longitude + eqn of equinoxes + sidereal DUT (radians)
*       (14)     local apparent sidereal time (radians)
*
*  Returned:
*     RAP    d      geocentric apparent right ascension
*     DAP    d      geocentric apparent declination
*
*  Notes:
*
*  1)  Only the first character of the TYPE argument is significant.
*      'R' or 'r' indicates that OBS1 and OBS2 are the observed right
*      ascension and declination;  'H' or 'h' indicates that they are
*      hour angle (west +ve) and declination;  anything else ('A' or
*      'a' is recommended) indicates that OBS1 and OBS2 are azimuth
*      (north zero, east 90 deg) and zenith distance.  (Zenith distance
*      is used rather than elevation in order to reflect the fact that
*      no allowance is made for depression of the horizon.)
*
*  2)  The accuracy of the result is limited by the corrections for
*      refraction.  Providing the meteorological parameters are
*      known accurately and there are no gross local effects, the
*      predicted apparent RA,Dec should be within about 0.1 arcsec
*      for a zenith distance of less than 70 degrees.  Even at a
*      topocentric zenith distance of 90 degrees, the accuracy in
*      elevation should be better than 1 arcmin;  useful results
*      are available for a further 3 degrees, beyond which the
*      sla_REFRO routine returns a fixed value of the refraction.
*      The complementary routines sla_AOP (or sla_AOPQK) and sla_OAP
*      (or sla_OAPQK) are self-consistent to better than 1 micro-
*      arcsecond all over the celestial sphere.
*
*  3)  It is advisable to take great care with units, as even
*      unlikely values of the input parameters are accepted and
*      processed in accordance with the models used.
*
*  5)  "Observed" Az,El means the position that would be seen by a
*      perfect theodolite located at the observer.  This is
*      related to the observed HA,Dec via the standard rotation, using
*      the geodetic latitude (corrected for polar motion), while the
*      observed HA and RA are related simply through the local
*      apparent ST.  "Observed" RA,Dec or HA,Dec thus means the
*      position that would be seen by a perfect equatorial located
*      at the observer and with its polar axis aligned to the
*      Earth's axis of rotation (n.b. not to the refracted pole).
*      By removing from the observed place the effects of
*      atmospheric refraction and diurnal aberration, the
*      geocentric apparent RA,Dec is obtained.
*
*  5)  Frequently, mean rather than apparent RA,Dec will be required,
*      in which case further transformations will be necessary.  The
*      sla_AMP etc routines will convert the apparent RA,Dec produced
*      by the present routine into an "FK5" (J2000) mean place, by
*      allowing for the Sun's gravitational lens effect, annual
*      aberration, nutation and precession.  Should "FK4" (1950)
*      coordinates be needed, the routines sla_FK524 etc will also
*      need to be applied.
*
*  6)  To convert to apparent RA,Dec the coordinates read from a
*      real telescope, corrections would have to be applied for
*      encoder zero points, gear and encoder errors, tube flexure,
*      the position of the rotator axis and the pointing axis
*      relative to it, non-perpendicularity between the mounting
*      axes, and finally for the tilt of the azimuth or polar axis
*      of the mounting (with appropriate corrections for mount
*      flexures).  Some telescopes would, of course, exhibit other
*      properties which would need to be accounted for at the
*      appropriate point in the sequence.
*
*  7)  The star-independent apparent-to-observed-place parameters
*      in AOPRMS may be computed by means of the sla_AOPPA routine.
*      If nothing has changed significantly except the time, the
*      sla_AOPPAT routine may be used to perform the requisite
*      partial recomputation of AOPRMS.
*
*  8) The azimuths etc used by the present routine are with respect
*     to the celestial pole.  Corrections from the terrestrial pole
*     can be computed using sla_POLMO.
*
*  Called:  sla_DCS2C, sla_DCC2S, sla_REFRO, sla_DRANRM
*
*  Last revision:   29 December 2004
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

      CHARACTER*(*) TYPE
      DOUBLE PRECISION OB1,OB2,AOPRMS(14),RAP,DAP

*  Breakpoint for fast/slow refraction algorithm:
*  ZD greater than arctan(4), (see sla_REFCO routine)
*  or vector Z less than cosine(arctan(Z)) = 1/sqrt(17)
      DOUBLE PRECISION ZBREAK
      PARAMETER (ZBREAK=0.242535625D0)

      CHARACTER C
      DOUBLE PRECISION C1,C2,SPHI,CPHI,ST,CE,XAEO,YAEO,ZAEO,V(3),
     :                 XMHDO,YMHDO,ZMHDO,AZ,SZ,ZDO,TZ,DREF,ZDT,
     :                 XAET,YAET,ZAET,XMHDA,YMHDA,ZMHDA,DIURAB,F,HMA

      DOUBLE PRECISION sla_DRANRM



*  Coordinate type
      C = TYPE(1:1)

*  Coordinates
      C1 = OB1
      C2 = OB2

*  Sin, cos of latitude
      SPHI = AOPRMS(2)
      CPHI = AOPRMS(3)

*  Local apparent sidereal time
      ST = AOPRMS(14)

*  Standardise coordinate type
      IF (C.EQ.'R'.OR.C.EQ.'r') THEN
         C = 'R'
      ELSE IF (C.EQ.'H'.OR.C.EQ.'h') THEN
         C = 'H'
      ELSE
         C = 'A'
      END IF

*  If Az,ZD convert to Cartesian (S=0,E=90)
      IF (C.EQ.'A') THEN
         CE = SIN(C2)
         XAEO = -COS(C1)*CE
         YAEO = SIN(C1)*CE
         ZAEO = COS(C2)
      ELSE

*     If RA,Dec convert to HA,Dec
         IF (C.EQ.'R') THEN
            C1 = ST-C1
         END IF

*     To Cartesian -HA,Dec
         CALL sla_DCS2C(-C1,C2,V)
         XMHDO = V(1)
         YMHDO = V(2)
         ZMHDO = V(3)

*     To Cartesian Az,El (S=0,E=90)
         XAEO = SPHI*XMHDO-CPHI*ZMHDO
         YAEO = YMHDO
         ZAEO = CPHI*XMHDO+SPHI*ZMHDO
      END IF

*  Azimuth (S=0,E=90)
      IF (XAEO.NE.0D0.OR.YAEO.NE.0D0) THEN
         AZ = ATAN2(YAEO,XAEO)
      ELSE
         AZ = 0D0
      END IF

*  Sine of observed ZD, and observed ZD
      SZ = SQRT(XAEO*XAEO+YAEO*YAEO)
      ZDO = ATAN2(SZ,ZAEO)

*
*  Refraction
*  ----------

*  Large zenith distance?
      IF (ZAEO.GE.ZBREAK) THEN

*     Fast algorithm using two constant model
         TZ = SZ/ZAEO
         DREF = (AOPRMS(11)+AOPRMS(12)*TZ*TZ)*TZ

      ELSE

*     Rigorous algorithm for large ZD
         CALL sla_REFRO(ZDO,AOPRMS(5),AOPRMS(6),AOPRMS(7),AOPRMS(8),
     :                  AOPRMS(9),AOPRMS(1),AOPRMS(10),1D-8,DREF)
      END IF

      ZDT = ZDO+DREF

*  To Cartesian Az,ZD
      CE = SIN(ZDT)
      XAET = COS(AZ)*CE
      YAET = SIN(AZ)*CE
      ZAET = COS(ZDT)

*  Cartesian Az,ZD to Cartesian -HA,Dec
      XMHDA = SPHI*XAET+CPHI*ZAET
      YMHDA = YAET
      ZMHDA = -CPHI*XAET+SPHI*ZAET

*  Diurnal aberration
      DIURAB = -AOPRMS(4)
      F = (1D0-DIURAB*YMHDA)
      V(1) = F*XMHDA
      V(2) = F*(YMHDA+DIURAB)
      V(3) = F*ZMHDA

*  To spherical -HA,Dec
      CALL sla_DCC2S(V,HMA,DAP)

*  Right Ascension
      RAP = sla_DRANRM(ST+HMA)

      END
