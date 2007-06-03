      SUBROUTINE sla_AOPPA ( DATE, DUT, ELONGM, PHIM, HM,
     :                       XP, YP, TDK, PMB, RH, WL, TLR, AOPRMS )
*+
*     - - - - - -
*      A O P P A
*     - - - - - -
*
*  Precompute apparent to observed place parameters required by
*  sla_AOPQK and sla_OAPQK.
*
*  Given:
*     DATE   d      UTC date/time (modified Julian Date, JD-2400000.5)
*     DUT    d      delta UT:  UT1-UTC (UTC seconds)
*     ELONGM d      mean longitude of the observer (radians, east +ve)
*     PHIM   d      mean geodetic latitude of the observer (radians)
*     HM     d      observer's height above sea level (metres)
*     XP     d      polar motion x-coordinate (radians)
*     YP     d      polar motion y-coordinate (radians)
*     TDK    d      local ambient temperature (K; std=273.15D0)
*     PMB    d      local atmospheric pressure (mb; std=1013.25D0)
*     RH     d      local relative humidity (in the range 0D0-1D0)
*     WL     d      effective wavelength (micron, e.g. 0.55D0)
*     TLR    d      tropospheric lapse rate (K/metre, e.g. 0.0065D0)
*
*  Returned:
*     AOPRMS d(14)  star-independent apparent-to-observed parameters:
*
*       (1)      geodetic latitude (radians)
*       (2,3)    sine and cosine of geodetic latitude
*       (4)      magnitude of diurnal aberration vector
*       (5)      height (HM)
*       (6)      ambient temperature (TDK)
*       (7)      pressure (PMB)
*       (8)      relative humidity (RH)
*       (9)      wavelength (WL)
*       (10)     lapse rate (TLR)
*       (11,12)  refraction constants A and B (radians)
*       (13)     longitude + eqn of equinoxes + sidereal DUT (radians)
*       (14)     local apparent sidereal time (radians)
*
*  Notes:
*
*   1)  It is advisable to take great care with units, as even
*       unlikely values of the input parameters are accepted and
*       processed in accordance with the models used.
*
*   2)  The DATE argument is UTC expressed as an MJD.  This is,
*       strictly speaking, improper, because of leap seconds.  However,
*       as long as the delta UT and the UTC are consistent there
*       are no difficulties, except during a leap second.  In this
*       case, the start of the 61st second of the final minute should
*       begin a new MJD day and the old pre-leap delta UT should
*       continue to be used.  As the 61st second completes, the MJD
*       should revert to the start of the day as, simultaneously,
*       the delta UTC changes by one second to its post-leap new value.
*
*   3)  The delta UT (UT1-UTC) is tabulated in IERS circulars and
*       elsewhere.  It increases by exactly one second at the end of
*       each UTC leap second, introduced in order to keep delta UT
*       within +/- 0.9 seconds.
*
*   4)  IMPORTANT -- TAKE CARE WITH THE LONGITUDE SIGN CONVENTION.
*       The longitude required by the present routine is east-positive,
*       in accordance with geographical convention (and right-handed).
*       In particular, note that the longitudes returned by the
*       sla_OBS routine are west-positive, following astronomical
*       usage, and must be reversed in sign before use in the present
*       routine.
*
*   5)  The polar coordinates XP,YP can be obtained from IERS
*       circulars and equivalent publications.  The maximum amplitude
*       is about 0.3 arcseconds.  If XP,YP values are unavailable,
*       use XP=YP=0D0.  See page B60 of the 1988 Astronomical Almanac
*       for a definition of the two angles.
*
*   6)  The height above sea level of the observing station, HM,
*       can be obtained from the Astronomical Almanac (Section J
*       in the 1988 edition), or via the routine sla_OBS.  If P,
*       the pressure in millibars, is available, an adequate
*       estimate of HM can be obtained from the expression
*
*             HM ~ -29.3D0*TSL*LOG(P/1013.25D0).
*
*       where TSL is the approximate sea-level air temperature in K
*       (see Astrophysical Quantities, C.W.Allen, 3rd edition,
*       section 52).  Similarly, if the pressure P is not known,
*       it can be estimated from the height of the observing
*       station, HM, as follows:
*
*             P ~ 1013.25D0*EXP(-HM/(29.3D0*TSL)).
*
*       Note, however, that the refraction is nearly proportional to the
*       pressure and that an accurate P value is important for precise
*       work.
*
*   7)  Repeated, computationally-expensive, calls to sla_AOPPA for
*       times that are very close together can be avoided by calling
*       sla_AOPPA just once and then using sla_AOPPAT for the subsequent
*       times.  Fresh calls to sla_AOPPA will be needed only when
*       changes in the precession have grown to unacceptable levels or
*       when anything affecting the refraction has changed.
*
*  Called:  sla_GEOC, sla_REFCO, sla_EQEQX, sla_AOPPAT
*
*  Last revision:   2 December 2005
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

      DOUBLE PRECISION DATE,DUT,ELONGM,PHIM,HM,XP,YP,TDK,PMB,
     :                 RH,WL,TLR,AOPRMS(14)

      DOUBLE PRECISION sla_EQEQX

*  2Pi
      DOUBLE PRECISION D2PI
      PARAMETER (D2PI=6.283185307179586476925287D0)

*  Seconds of time to radians
      DOUBLE PRECISION S2R
      PARAMETER (S2R=7.272205216643039903848712D-5)

*  Speed of light (AU per day)
      DOUBLE PRECISION C
      PARAMETER (C=173.14463331D0)

*  Ratio between solar and sidereal time
      DOUBLE PRECISION SOLSID
      PARAMETER (SOLSID=1.00273790935D0)

      DOUBLE PRECISION CPHIM,XT,YT,ZT,XC,YC,ZC,ELONG,PHI,UAU,VAU



*  Observer's location corrected for polar motion
      CPHIM = COS(PHIM)
      XT = COS(ELONGM)*CPHIM
      YT = SIN(ELONGM)*CPHIM
      ZT = SIN(PHIM)
      XC = XT-XP*ZT
      YC = YT+YP*ZT
      ZC = XP*XT-YP*YT+ZT
      IF (XC.EQ.0D0.AND.YC.EQ.0D0) THEN
         ELONG = 0D0
      ELSE
         ELONG = ATAN2(YC,XC)
      END IF
      PHI = ATAN2(ZC,SQRT(XC*XC+YC*YC))
      AOPRMS(1) = PHI
      AOPRMS(2) = SIN(PHI)
      AOPRMS(3) = COS(PHI)

*  Magnitude of the diurnal aberration vector
      CALL sla_GEOC(PHI,HM,UAU,VAU)
      AOPRMS(4) = D2PI*UAU*SOLSID/C

*  Copy the refraction parameters and compute the A & B constants
      AOPRMS(5) = HM
      AOPRMS(6) = TDK
      AOPRMS(7) = PMB
      AOPRMS(8) = RH
      AOPRMS(9) = WL
      AOPRMS(10) = TLR
      CALL sla_REFCO(HM,TDK,PMB,RH,WL,PHI,TLR,1D-10,
     :               AOPRMS(11),AOPRMS(12))

*  Longitude + equation of the equinoxes + sidereal equivalent of DUT
*  (ignoring change in equation of the equinoxes between UTC and TDB)
      AOPRMS(13) = ELONG+sla_EQEQX(DATE)+DUT*SOLSID*S2R

*  Sidereal time
      CALL sla_AOPPAT(DATE,AOPRMS)

      END
