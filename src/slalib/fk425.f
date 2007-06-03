      SUBROUTINE sla_FK425 (R1950,D1950,DR1950,DD1950,P1950,V1950,
     :                      R2000,D2000,DR2000,DD2000,P2000,V2000)
*+
*     - - - - - -
*      F K 4 2 5
*     - - - - - -
*
*  Convert B1950.0 FK4 star data to J2000.0 FK5 (double precision)
*
*  This routine converts stars from the old, Bessel-Newcomb, FK4
*  system to the new, IAU 1976, FK5, Fricke system.  The precepts
*  of Smith et al (Ref 1) are followed, using the implementation
*  by Yallop et al (Ref 2) of a matrix method due to Standish.
*  Kinoshita's development of Andoyer's post-Newcomb precession is
*  used.  The numerical constants from Seidelmann et al (Ref 3) are
*  used canonically.
*
*  Given:  (all B1950.0,FK4)
*     R1950,D1950     dp    B1950.0 RA,Dec (rad)
*     DR1950,DD1950   dp    B1950.0 proper motions (rad/trop.yr)
*     P1950           dp    parallax (arcsec)
*     V1950           dp    radial velocity (km/s, +ve = moving away)
*
*  Returned:  (all J2000.0,FK5)
*     R2000,D2000     dp    J2000.0 RA,Dec (rad)
*     DR2000,DD2000   dp    J2000.0 proper motions (rad/Jul.yr)
*     P2000           dp    parallax (arcsec)
*     V2000           dp    radial velocity (km/s, +ve = moving away)
*
*  Notes:
*
*  1)  The proper motions in RA are dRA/dt rather than
*      cos(Dec)*dRA/dt, and are per year rather than per century.
*
*  2)  Conversion from Besselian epoch 1950.0 to Julian epoch
*      2000.0 only is provided for.  Conversions involving other
*      epochs will require use of the appropriate precession,
*      proper motion, and E-terms routines before and/or
*      after FK425 is called.
*
*  3)  In the FK4 catalogue the proper motions of stars within
*      10 degrees of the poles do not embody the differential
*      E-term effect and should, strictly speaking, be handled
*      in a different manner from stars outside these regions.
*      However, given the general lack of homogeneity of the star
*      data available for routine astrometry, the difficulties of
*      handling positions that may have been determined from
*      astrometric fields spanning the polar and non-polar regions,
*      the likelihood that the differential E-terms effect was not
*      taken into account when allowing for proper motion in past
*      astrometry, and the undesirability of a discontinuity in
*      the algorithm, the decision has been made in this routine to
*      include the effect of differential E-terms on the proper
*      motions for all stars, whether polar or not.  At epoch 2000,
*      and measuring on the sky rather than in terms of dRA, the
*      errors resulting from this simplification are less than
*      1 milliarcsecond in position and 1 milliarcsecond per
*      century in proper motion.
*
*  References:
*
*     1  Smith, C.A. et al, 1989.  "The transformation of astrometric
*        catalog systems to the equinox J2000.0".  Astron.J. 97, 265.
*
*     2  Yallop, B.D. et al, 1989.  "Transformation of mean star places
*        from FK4 B1950.0 to FK5 J2000.0 using matrices in 6-space".
*        Astron.J. 97, 274.
*
*     3  Seidelmann, P.K. (ed), 1992.  "Explanatory Supplement to
*        the Astronomical Almanac", ISBN 0-935702-68-7.
*
*  P.T.Wallace   Starlink   19 December 1993
*
*  Copyright (C) 1995 Rutherford Appleton Laboratory
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

      DOUBLE PRECISION R1950,D1950,DR1950,DD1950,P1950,V1950,
     :                 R2000,D2000,DR2000,DD2000,P2000,V2000


*  Miscellaneous
      DOUBLE PRECISION R,D,UR,UD,PX,RV,SR,CR,SD,CD,W,WD
      DOUBLE PRECISION X,Y,Z,XD,YD,ZD
      DOUBLE PRECISION RXYSQ,RXYZSQ,RXY,RXYZ,SPXY,SPXYZ
      INTEGER I,J

*  Star position and velocity vectors
      DOUBLE PRECISION R0(3),RD0(3)

*  Combined position and velocity vectors
      DOUBLE PRECISION V1(6),V2(6)

*  2Pi
      DOUBLE PRECISION D2PI
      PARAMETER (D2PI=6.283185307179586476925287D0)

*  Radians per year to arcsec per century
      DOUBLE PRECISION PMF
      PARAMETER (PMF=100D0*60D0*60D0*360D0/D2PI)

*  Small number to avoid arithmetic problems
      DOUBLE PRECISION TINY
      PARAMETER (TINY=1D-30)


*
*  CANONICAL CONSTANTS  (see references)
*

*  Km per sec to AU per tropical century
*  = 86400 * 36524.2198782 / 149597870
      DOUBLE PRECISION VF
      PARAMETER (VF=21.095D0)

*  Constant vector and matrix (by columns)
      DOUBLE PRECISION A(3),AD(3),EM(6,6)
      DATA A,AD/ -1.62557D-6,  -0.31919D-6, -0.13843D-6,
     :           +1.245D-3,    -1.580D-3,   -0.659D-3/

      DATA (EM(I,1),I=1,6) / +0.9999256782D0,
     :                       +0.0111820610D0,
     :                       +0.0048579479D0,
     :                       -0.000551D0,
     :                       +0.238514D0,
     :                       -0.435623D0 /

      DATA (EM(I,2),I=1,6) / -0.0111820611D0,
     :                       +0.9999374784D0,
     :                       -0.0000271474D0,
     :                       -0.238565D0,
     :                       -0.002667D0,
     :                       +0.012254D0 /

      DATA (EM(I,3),I=1,6) / -0.0048579477D0,
     :                       -0.0000271765D0,
     :                       +0.9999881997D0,
     :                       +0.435739D0,
     :                       -0.008541D0,
     :                       +0.002117D0 /

      DATA (EM(I,4),I=1,6) / +0.00000242395018D0,
     :                       +0.00000002710663D0,
     :                       +0.00000001177656D0,
     :                       +0.99994704D0,
     :                       +0.01118251D0,
     :                       +0.00485767D0 /

      DATA (EM(I,5),I=1,6) / -0.00000002710663D0,
     :                       +0.00000242397878D0,
     :                       -0.00000000006582D0,
     :                       -0.01118251D0,
     :                       +0.99995883D0,
     :                       -0.00002714D0 /

      DATA (EM(I,6),I=1,6) / -0.00000001177656D0,
     :                       -0.00000000006587D0,
     :                       +0.00000242410173D0,
     :                       -0.00485767D0,
     :                       -0.00002718D0,
     :                       +1.00000956D0 /



*  Pick up B1950 data (units radians and arcsec/TC)
      R=R1950
      D=D1950
      UR=DR1950*PMF
      UD=DD1950*PMF
      PX=P1950
      RV=V1950

*  Spherical to Cartesian
      SR=SIN(R)
      CR=COS(R)
      SD=SIN(D)
      CD=COS(D)

      R0(1)=CR*CD
      R0(2)=SR*CD
      R0(3)=   SD

      W=VF*RV*PX

      RD0(1)=-SR*CD*UR-CR*SD*UD+W*R0(1)
      RD0(2)= CR*CD*UR-SR*SD*UD+W*R0(2)
      RD0(3)=             CD*UD+W*R0(3)

*  Allow for e-terms and express as position+velocity 6-vector
      W=R0(1)*A(1)+R0(2)*A(2)+R0(3)*A(3)
      WD=R0(1)*AD(1)+R0(2)*AD(2)+R0(3)*AD(3)
      DO I=1,3
         V1(I)=R0(I)-A(I)+W*R0(I)
         V1(I+3)=RD0(I)-AD(I)+WD*R0(I)
      END DO

*  Convert position+velocity vector to Fricke system
      DO I=1,6
         W=0D0
         DO J=1,6
            W=W+EM(I,J)*V1(J)
         END DO
         V2(I)=W
      END DO

*  Revert to spherical coordinates
      X=V2(1)
      Y=V2(2)
      Z=V2(3)
      XD=V2(4)
      YD=V2(5)
      ZD=V2(6)

      RXYSQ=X*X+Y*Y
      RXYZSQ=RXYSQ+Z*Z
      RXY=SQRT(RXYSQ)
      RXYZ=SQRT(RXYZSQ)

      SPXY=X*XD+Y*YD
      SPXYZ=SPXY+Z*ZD

      IF (X.EQ.0D0.AND.Y.EQ.0D0) THEN
         R=0D0
      ELSE
         R=ATAN2(Y,X)
         IF (R.LT.0.0D0) R=R+D2PI
      END IF
      D=ATAN2(Z,RXY)

      IF (RXY.GT.TINY) THEN
         UR=(X*YD-Y*XD)/RXYSQ
         UD=(ZD*RXYSQ-Z*SPXY)/(RXYZSQ*RXY)
      END IF

      IF (PX.GT.TINY) THEN
         RV=SPXYZ/(PX*RXYZ*VF)
         PX=PX/RXYZ
      END IF

*  Return results
      R2000=R
      D2000=D
      DR2000=UR/PMF
      DD2000=UD/PMF
      V2000=RV
      P2000=PX

      END
