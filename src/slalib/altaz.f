      SUBROUTINE sla_ALTAZ (HA, DEC, PHI,
     :                      AZ, AZD, AZDD, EL, ELD, ELDD, PA, PAD, PADD)
*+
*     - - - - - -
*      A L T A Z
*     - - - - - -
*
*  Positions, velocities and accelerations for an altazimuth
*  telescope mount.
*
*  (double precision)
*
*  Given:
*     HA      d     hour angle
*     DEC     d     declination
*     PHI     d     observatory latitude
*
*  Returned:
*     AZ      d     azimuth
*     AZD     d        "    velocity
*     AZDD    d        "    acceleration
*     EL      d     elevation
*     ELD     d         "     velocity
*     ELDD    d         "     acceleration
*     PA      d     parallactic angle
*     PAD     d         "      "   velocity
*     PADD    d         "      "   acceleration
*
*  Notes:
*
*  1)  Natural units are used throughout.  HA, DEC, PHI, AZ, EL
*      and ZD are in radians.  The velocities and accelerations
*      assume constant declination and constant rate of change of
*      hour angle (as for tracking a star);  the units of AZD, ELD
*      and PAD are radians per radian of HA, while the units of AZDD,
*      ELDD and PADD are radians per radian of HA squared.  To
*      convert into practical degree- and second-based units:
*
*        angles * 360/2pi -> degrees
*        velocities * (2pi/86400)*(360/2pi) -> degree/sec
*        accelerations * ((2pi/86400)**2)*(360/2pi) -> degree/sec/sec
*
*      Note that the seconds here are sidereal rather than SI.  One
*      sidereal second is about 0.99727 SI seconds.
*
*      The velocity and acceleration factors assume the sidereal
*      tracking case.  Their respective numerical values are (exactly)
*      1/240 and (approximately) 1/3300236.9.
*
*  2)  Azimuth is returned in the range 0-2pi;  north is zero,
*      and east is +pi/2.  Elevation and parallactic angle are
*      returned in the range +/-pi.  Parallactic angle is +ve for
*      a star west of the meridian and is the angle NP-star-zenith.
*
*  3)  The latitude is geodetic as opposed to geocentric.  The
*      hour angle and declination are topocentric.  Refraction and
*      deficiencies in the telescope mounting are ignored.  The
*      purpose of the routine is to give the general form of the
*      quantities.  The details of a real telescope could profoundly
*      change the results, especially close to the zenith.
*
*  4)  No range checking of arguments is carried out.
*
*  5)  In applications which involve many such calculations, rather
*      than calling the present routine it will be more efficient to
*      use inline code, having previously computed fixed terms such
*      as sine and cosine of latitude, and (for tracking a star)
*      sine and cosine of declination.
*
*  This revision:  29 October 2004
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

      DOUBLE PRECISION HA,DEC,PHI,AZ,AZD,AZDD,EL,ELD,ELDD,PA,PAD,PADD

      DOUBLE PRECISION DPI,D2PI,TINY
      PARAMETER (DPI=3.1415926535897932384626433832795D0,
     :           D2PI=6.283185307179586476925286766559D0,
     :           TINY=1D-30)

      DOUBLE PRECISION SH,CH,SD,CD,SP,CP,CHCD,SDCP,X,Y,Z,RSQ,R,A,E,C,S,
     :                 Q,QD,AD,ED,EDR,ADD,EDD,QDD


*  Useful functions
      SH=SIN(HA)
      CH=COS(HA)
      SD=SIN(DEC)
      CD=COS(DEC)
      SP=SIN(PHI)
      CP=COS(PHI)
      CHCD=CH*CD
      SDCP=SD*CP
      X=-CHCD*SP+SDCP
      Y=-SH*CD
      Z=CHCD*CP+SD*SP
      RSQ=X*X+Y*Y
      R=SQRT(RSQ)

*  Azimuth and elevation
      IF (RSQ.EQ.0D0) THEN
         A=0D0
      ELSE
         A=ATAN2(Y,X)
      END IF
      IF (A.LT.0D0) A=A+D2PI
      E=ATAN2(Z,R)

*  Parallactic angle
      C=CD*SP-CH*SDCP
      S=SH*CP
      IF (C*C+S*S.GT.0) THEN
         Q=ATAN2(S,C)
      ELSE
         Q=DPI-HA
      END IF

*  Velocities and accelerations (clamped at zenith/nadir)
      IF (RSQ.LT.TINY) THEN
         RSQ=TINY
         R=SQRT(RSQ)
      END IF
      QD=-X*CP/RSQ
      AD=SP+Z*QD
      ED=CP*Y/R
      EDR=ED/R
      ADD=EDR*(Z*SP+(2D0-RSQ)*QD)
      EDD=-R*QD*AD
      QDD=EDR*(SP+2D0*Z*QD)

*  Results
      AZ=A
      AZD=AD
      AZDD=ADD
      EL=E
      ELD=ED
      ELDD=EDD
      PA=Q
      PAD=QD
      PADD=QDD

      END
