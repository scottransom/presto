      SUBROUTINE sla_AOPQK (RAP, DAP, AOPRMS, AOB, ZOB, HOB, DOB, ROB)
*+
*     - - - - - -
*      A O P Q K
*     - - - - - -
*
*  Quick apparent to observed place (but see note 8, below, for
*  remarks about speed).
*
*  Given:
*     RAP    d      geocentric apparent right ascension
*     DAP    d      geocentric apparent declination
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
*     AOB    d      observed azimuth (radians: N=0,E=90)
*     ZOB    d      observed zenith distance (radians)
*     HOB    d      observed Hour Angle (radians)
*     DOB    d      observed Declination (radians)
*     ROB    d      observed Right Ascension (radians)
*
*  Notes:
*
*   1)  This routine returns zenith distance rather than elevation
*       in order to reflect the fact that no allowance is made for
*       depression of the horizon.
*
*   2)  The accuracy of the result is limited by the corrections for
*       refraction.  Providing the meteorological parameters are
*       known accurately and there are no gross local effects, the
*       observed RA,Dec predicted by this routine should be within
*       about 0.1 arcsec for a zenith distance of less than 70 degrees.
*       Even at a topocentric zenith distance of 90 degrees, the
*       accuracy in elevation should be better than 1 arcmin;  useful
*       results are available for a further 3 degrees, beyond which
*       the sla_REFRO routine returns a fixed value of the refraction.
*       The complementary routines sla_AOP (or sla_AOPQK) and sla_OaAP
*       (or sla_OAPQK) are self-consistent to better than 1 micro-
*       arcsecond all over the celestial sphere.
*
*   3)  It is advisable to take great care with units, as even
*       unlikely values of the input parameters are accepted and
*       processed in accordance with the models used.
*
*   4)  "Apparent" place means the geocentric apparent right ascension
*       and declination, which is obtained from a catalogue mean place
*       by allowing for space motion, parallax, precession, nutation,
*       annual aberration, and the Sun's gravitational lens effect.  For
*       star positions in the FK5 system (i.e. J2000), these effects can
*       be applied by means of the sla_MAP etc routines.  Starting from
*       other mean place systems, additional transformations will be
*       needed;  for example, FK4 (i.e. B1950) mean places would first
*       have to be converted to FK5, which can be done with the
*       sla_FK425 etc routines.
*
*   5)  "Observed" Az,El means the position that would be seen by a
*       perfect theodolite located at the observer.  This is obtained
*       from the geocentric apparent RA,Dec by allowing for Earth
*       orientation and diurnal aberration, rotating from equator
*       to horizon coordinates, and then adjusting for refraction.
*       The HA,Dec is obtained by rotating back into equatorial
*       coordinates, using the geodetic latitude corrected for polar
*       motion, and is the position that would be seen by a perfect
*       equatorial located at the observer and with its polar axis
*       aligned to the Earth's axis of rotation (n.b. not to the
*       refracted pole).  Finally, the RA is obtained by subtracting
*       the HA from the local apparent ST.
*
*   6)  To predict the required setting of a real telescope, the
*       observed place produced by this routine would have to be
*       adjusted for the tilt of the azimuth or polar axis of the
*       mounting (with appropriate corrections for mount flexures),
*       for non-perpendicularity between the mounting axes, for the
*       position of the rotator axis and the pointing axis relative
*       to it, for tube flexure, for gear and encoder errors, and
*       finally for encoder zero points.  Some telescopes would, of
*       course, exhibit other properties which would need to be
*       accounted for at the appropriate point in the sequence.
*
*   7)  The star-independent apparent-to-observed-place parameters
*       in AOPRMS may be computed by means of the sla_AOPPA routine.
*       If nothing has changed significantly except the time, the
*       sla_AOPPAT routine may be used to perform the requisite
*       partial recomputation of AOPRMS.
*
*   8)  At zenith distances beyond about 76 degrees, the need for
*       special care with the corrections for refraction causes a
*       marked increase in execution time.  Moreover, the effect
*       gets worse with increasing zenith distance.  Adroit
*       programming in the calling application may allow the
*       problem to be reduced.  Prepare an alternative AOPRMS array,
*       computed for zero air-pressure;  this will disable the
*       refraction corrections and cause rapid execution.  Using
*       this AOPRMS array, a preliminary call to the present routine
*       will, depending on the application, produce a rough position
*       which may be enough to establish whether the full, slow
*       calculation (using the real AOPRMS array) is worthwhile.
*       For example, there would be no need for the full calculation
*       if the preliminary call had already established that the
*       source was well below the elevation limits for a particular
*       telescope.
*
*  9)   The azimuths etc produced by the present routine are with
*       respect to the celestial pole.  Corrections to the terrestrial
*       pole can be computed using sla_POLMO.
*
*  Called:  sla_DCS2C, sla_REFZ, sla_REFRO, sla_DCC2S, sla_DRANRM
*
*  P.T.Wallace   Starlink   24 October 2003
*
*  Copyright (C) 2003 Rutherford Appleton Laboratory
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

      DOUBLE PRECISION RAP,DAP,AOPRMS(14),AOB,ZOB,HOB,DOB,ROB

*  Breakpoint for fast/slow refraction algorithm:
*  ZD greater than arctan(4), (see sla_REFCO routine)
*  or vector Z less than cosine(arctan(Z)) = 1/sqrt(17)
      DOUBLE PRECISION ZBREAK
      PARAMETER (ZBREAK=0.242535625D0)

      INTEGER I

      DOUBLE PRECISION SPHI,CPHI,ST,V(3),XHD,YHD,ZHD,DIURAB,F,
     :                 XHDT,YHDT,ZHDT,XAET,YAET,ZAET,AZOBS,
     :                 ZDT,REFA,REFB,ZDOBS,DZD,DREF,CE,
     :                 XAEO,YAEO,ZAEO,HMOBS,DCOBS,RAOBS

      DOUBLE PRECISION sla_DRANRM



*  Sin, cos of latitude
      SPHI = AOPRMS(2)
      CPHI = AOPRMS(3)

*  Local apparent sidereal time
      ST = AOPRMS(14)

*  Apparent RA,Dec to Cartesian -HA,Dec
      CALL sla_DCS2C(RAP-ST,DAP,V)
      XHD = V(1)
      YHD = V(2)
      ZHD = V(3)

*  Diurnal aberration
      DIURAB = AOPRMS(4)
      F = (1D0-DIURAB*YHD)
      XHDT = F*XHD
      YHDT = F*(YHD+DIURAB)
      ZHDT = F*ZHD

*  Cartesian -HA,Dec to Cartesian Az,El (S=0,E=90)
      XAET = SPHI*XHDT-CPHI*ZHDT
      YAET = YHDT
      ZAET = CPHI*XHDT+SPHI*ZHDT

*  Azimuth (N=0,E=90)
      IF (XAET.EQ.0D0.AND.YAET.EQ.0D0) THEN
         AZOBS = 0D0
      ELSE
         AZOBS = ATAN2(YAET,-XAET)
      END IF

*  Topocentric zenith distance
      ZDT = ATAN2(SQRT(XAET*XAET+YAET*YAET),ZAET)

*
*  Refraction
*  ----------

*  Fast algorithm using two constant model
      REFA = AOPRMS(11)
      REFB = AOPRMS(12)
      CALL sla_REFZ(ZDT,REFA,REFB,ZDOBS)

*  Large zenith distance?
      IF (COS(ZDOBS).LT.ZBREAK) THEN

*     Yes: use rigorous algorithm

*     Initialize loop (maximum of 10 iterations)
         I = 1
         DZD = 1D1
         DO WHILE (ABS(DZD).GT.1D-10.AND.I.LE.10)

*        Compute refraction using current estimate of observed ZD
            CALL sla_REFRO(ZDOBS,AOPRMS(5),AOPRMS(6),AOPRMS(7),
     :                     AOPRMS(8),AOPRMS(9),AOPRMS(1),
     :                     AOPRMS(10),1D-8,DREF)

*        Remaining discrepancy
            DZD = ZDOBS+DREF-ZDT

*        Update the estimate
            ZDOBS = ZDOBS-DZD

*        Increment the iteration counter
            I = I+1
         END DO
      END IF

*  To Cartesian Az/ZD
      CE = SIN(ZDOBS)
      XAEO = -COS(AZOBS)*CE
      YAEO = SIN(AZOBS)*CE
      ZAEO = COS(ZDOBS)

*  Cartesian Az/ZD to Cartesian -HA,Dec
      V(1) = SPHI*XAEO+CPHI*ZAEO
      V(2) = YAEO
      V(3) = -CPHI*XAEO+SPHI*ZAEO

*  To spherical -HA,Dec
      CALL sla_DCC2S(V,HMOBS,DCOBS)

*  Right Ascension
      RAOBS = sla_DRANRM(ST+HMOBS)

*  Return the results
      AOB = AZOBS
      ZOB = ZDOBS
      HOB = -HMOBS
      DOB = DCOBS
      ROB = RAOBS

      END
