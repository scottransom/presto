      SUBROUTINE sla_OAP ( TYPE, OB1, OB2, DATE, DUT, ELONGM, PHIM,
     :                     HM, XP, YP, TDK, PMB, RH, WL, TLR,
     :                     RAP, DAP )
*+
*     - - - -
*      O A P
*     - - - -
*
*  Observed to apparent place.
*
*  Given:
*     TYPE   c*(*)  type of coordinates - 'R', 'H' or 'A' (see below)
*     OB1    d      observed Az, HA or RA (radians; Az is N=0,E=90)
*     OB2    d      observed ZD or Dec (radians)
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
*      (north zero, east 90 deg) and zenith distance.  (Zenith
*      distance is used rather than elevation in order to reflect the
*      fact that no allowance is made for depression of the horizon.)
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
*  4)  "Observed" Az,El means the position that would be seen by a
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
*  7)  This routine takes time to execute, due mainly to the rigorous
*      integration used to evaluate the refraction.  For processing
*      multiple stars for one location and time, call sla_AOPPA once
*      followed by one call per star to sla_OAPQK.  Where a range of
*      times within a limited period of a few hours is involved, and the
*      highest precision is not required, call sla_AOPPA once, followed
*      by a call to sla_AOPPAT each time the time changes, followed by
*      one call per star to sla_OAPQK.
*
*  8)  The DATE argument is UTC expressed as an MJD.  This is, strictly
*      speaking, wrong, because of leap seconds.  However, as long as
*      the delta UT and the UTC are consistent there are no
*      difficulties, except during a leap second.  In this case, the
*      start of the 61st second of the final minute should begin a new
*      MJD day and the old pre-leap delta UT should continue to be used.
*      As the 61st second completes, the MJD should revert to the start
*      of the day as, simultaneously, the delta UTC changes by one
*      second to its post-leap new value.
*
*  9)  The delta UT (UT1-UTC) is tabulated in IERS circulars and
*      elsewhere.  It increases by exactly one second at the end of
*      each UTC leap second, introduced in order to keep delta UT
*      within +/- 0.9 seconds.
*
*  10) IMPORTANT -- TAKE CARE WITH THE LONGITUDE SIGN CONVENTION.
*      The longitude required by the present routine is east-positive,
*      in accordance with geographical convention (and right-handed).
*      In particular, note that the longitudes returned by the
*      sla_OBS routine are west-positive, following astronomical
*      usage, and must be reversed in sign before use in the present
*      routine.
*
*  11) The polar coordinates XP,YP can be obtained from IERS
*      circulars and equivalent publications.  The maximum amplitude
*      is about 0.3 arcseconds.  If XP,YP values are unavailable,
*      use XP=YP=0D0.  See page B60 of the 1988 Astronomical Almanac
*      for a definition of the two angles.
*
*  12) The height above sea level of the observing station, HM,
*      can be obtained from the Astronomical Almanac (Section J
*      in the 1988 edition), or via the routine sla_OBS.  If P,
*      the pressure in millibars, is available, an adequate
*      estimate of HM can be obtained from the expression
*
*             HM ~ -29.3D0*TSL*LOG(P/1013.25D0).
*
*      where TSL is the approximate sea-level air temperature in K
*      (see Astrophysical Quantities, C.W.Allen, 3rd edition,
*      section 52).  Similarly, if the pressure P is not known,
*      it can be estimated from the height of the observing
*      station, HM, as follows:
*
*             P ~ 1013.25D0*EXP(-HM/(29.3D0*TSL)).
*
*      Note, however, that the refraction is nearly proportional to the
*      pressure and that an accurate P value is important for precise
*      work.
*
*  13) The azimuths etc. used by the present routine are with respect
*      to the celestial pole.  Corrections from the terrestrial pole
*      can be computed using sla_POLMO.
*
*  Called:  sla_AOPPA, sla_OAPQK
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

      CHARACTER*(*) TYPE
      DOUBLE PRECISION OB1,OB2,DATE,DUT,ELONGM,PHIM,HM,
     :                 XP,YP,TDK,PMB,RH,WL,TLR,RAP,DAP

      DOUBLE PRECISION AOPRMS(14)


      CALL sla_AOPPA(DATE,DUT,ELONGM,PHIM,HM,XP,YP,TDK,PMB,RH,WL,TLR,
     :               AOPRMS)
      CALL sla_OAPQK(TYPE,OB1,OB2,AOPRMS,RAP,DAP)

      END
