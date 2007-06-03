      SUBROUTINE sla_POLMO ( ELONGM, PHIM, XP, YP, ELONG, PHI, DAZ )
*+
*     - - - - - -
*      P O L M O
*     - - - - - -
*
*  Polar motion:  correct site longitude and latitude for polar
*  motion and calculate azimuth difference between celestial and
*  terrestrial poles.
*
*  Given:
*     ELONGM   d      mean longitude of the observer (radians, east +ve)
*     PHIM     d      mean geodetic latitude of the observer (radians)
*     XP       d      polar motion x-coordinate (radians)
*     YP       d      polar motion y-coordinate (radians)
*
*  Returned:
*     ELONG    d      true longitude of the observer (radians, east +ve)
*     PHI      d      true geodetic latitude of the observer (radians)
*     DAZ      d      azimuth correction (terrestrial-celestial, radians)
*
*  Notes:
*
*   1)  "Mean" longitude and latitude are the (fixed) values for the
*       site's location with respect to the IERS terrestrial reference
*       frame;  the latitude is geodetic.  TAKE CARE WITH THE LONGITUDE
*       SIGN CONVENTION.  The longitudes used by the present routine
*       are east-positive, in accordance with geographical convention
*       (and right-handed).  In particular, note that the longitudes
*       returned by the sla_OBS routine are west-positive, following
*       astronomical usage, and must be reversed in sign before use in
*       the present routine.
*
*   2)  XP and YP are the (changing) coordinates of the Celestial
*       Ephemeris Pole with respect to the IERS Reference Pole.
*       XP is positive along the meridian at longitude 0 degrees,
*       and YP is positive along the meridian at longitude
*       270 degrees (i.e. 90 degrees west).  Values for XP,YP can
*       be obtained from IERS circulars and equivalent publications;
*       the maximum amplitude observed so far is about 0.3 arcseconds.
*
*   3)  "True" longitude and latitude are the (moving) values for
*       the site's location with respect to the celestial ephemeris
*       pole and the meridian which corresponds to the Greenwich
*       apparent sidereal time.  The true longitude and latitude
*       link the terrestrial coordinates with the standard celestial
*       models (for precession, nutation, sidereal time etc).
*
*   4)  The azimuths produced by sla_AOP and sla_AOPQK are with
*       respect to due north as defined by the Celestial Ephemeris
*       Pole, and can therefore be called "celestial azimuths".
*       However, a telescope fixed to the Earth measures azimuth
*       essentially with respect to due north as defined by the
*       IERS Reference Pole, and can therefore be called "terrestrial
*       azimuth".  Uncorrected, this would manifest itself as a
*       changing "azimuth zero-point error".  The value DAZ is the
*       correction to be added to a celestial azimuth to produce
*       a terrestrial azimuth.
*
*   5)  The present routine is rigorous.  For most practical
*       purposes, the following simplified formulae provide an
*       adequate approximation:
*
*       ELONG = ELONGM+XP*COS(ELONGM)-YP*SIN(ELONGM)
*       PHI   = PHIM+(XP*SIN(ELONGM)+YP*COS(ELONGM))*TAN(PHIM)
*       DAZ   = -SQRT(XP*XP+YP*YP)*COS(ELONGM-ATAN2(XP,YP))/COS(PHIM)
*
*       An alternative formulation for DAZ is:
*
*       X = COS(ELONGM)*COS(PHIM)
*       Y = SIN(ELONGM)*COS(PHIM)
*       DAZ = ATAN2(-X*YP-Y*XP,X*X+Y*Y)
*
*   Reference:  Seidelmann, P.K. (ed), 1992.  "Explanatory Supplement
*               to the Astronomical Almanac", ISBN 0-935702-68-7,
*               sections 3.27, 4.25, 4.52.
*
*  P.T.Wallace   Starlink   30 November 2000
*
*  Copyright (C) 2000 Rutherford Appleton Laboratory
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

      DOUBLE PRECISION ELONGM,PHIM,XP,YP,ELONG,PHI,DAZ

      DOUBLE PRECISION SEL,CEL,SPH,CPH,XM,YM,ZM,XNM,YNM,ZNM,
     :                 SXP,CXP,SYP,CYP,ZW,XT,YT,ZT,XNT,YNT



*  Site mean longitude and mean geodetic latitude as a Cartesian vector
      SEL=SIN(ELONGM)
      CEL=COS(ELONGM)
      SPH=SIN(PHIM)
      CPH=COS(PHIM)

      XM=CEL*CPH
      YM=SEL*CPH
      ZM=SPH

*  Rotate site vector by polar motion, Y-component then X-component
      SXP=SIN(XP)
      CXP=COS(XP)
      SYP=SIN(YP)
      CYP=COS(YP)

      ZW=(-YM*SYP+ZM*CYP)

      XT=XM*CXP-ZW*SXP
      YT=YM*CYP+ZM*SYP
      ZT=XM*SXP+ZW*CXP

*  Rotate also the geocentric direction of the terrestrial pole (0,0,1)
      XNM=-SXP*CYP
      YNM=SYP
      ZNM=CXP*CYP

      CPH=SQRT(XT*XT+YT*YT)
      IF (CPH.EQ.0D0) XT=1D0
      SEL=YT/CPH
      CEL=XT/CPH

*  Return true longitude and true geodetic latitude of site
      IF (XT.NE.0D0.OR.YT.NE.0D0) THEN
         ELONG=ATAN2(YT,XT)
      ELSE
         ELONG=0D0
      END IF
      PHI=ATAN2(ZT,CPH)

*  Return current azimuth of terrestrial pole seen from site position
      XNT=(XNM*CEL+YNM*SEL)*ZT-ZNM*CPH
      YNT=-XNM*SEL+YNM*CEL
      IF (XNT.NE.0D0.OR.YNT.NE.0D0) THEN
         DAZ=ATAN2(-YNT,-XNT)
      ELSE
         DAZ=0D0
      END IF

      END
