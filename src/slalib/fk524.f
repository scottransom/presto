      SUBROUTINE sla_FK524 (R2000,D2000,DR2000,DD2000,P2000,V2000,
     :                      R1950,D1950,DR1950,DD1950,P1950,V1950)
*+
*     - - - - - -
*      F K 5 2 4
*     - - - - - -
*
*  Convert J2000.0 FK5 star data to B1950.0 FK4 (double precision)
*
*  This routine converts stars from the new, IAU 1976, FK5, Fricke
*  system, to the old, Bessel-Newcomb, FK4 system.  The precepts
*  of Smith et al (Ref 1) are followed, using the implementation
*  by Yallop et al (Ref 2) of a matrix method due to Standish.
*  Kinoshita's development of Andoyer's post-Newcomb precession is
*  used.  The numerical constants from Seidelmann et al (Ref 3) are
*  used canonically.
*
*  Given:  (all J2000.0,FK5)
*     R2000,D2000     dp    J2000.0 RA,Dec (rad)
*     DR2000,DD2000   dp    J2000.0 proper motions (rad/Jul.yr)
*     P2000           dp    parallax (arcsec)
*     V2000           dp    radial velocity (km/s, +ve = moving away)
*
*  Returned:  (all B1950.0,FK4)
*     R1950,D1950     dp    B1950.0 RA,Dec (rad)
*     DR1950,DD1950   dp    B1950.0 proper motions (rad/trop.yr)
*     P1950           dp    parallax (arcsec)
*     V1950           dp    radial velocity (km/s, +ve = moving away)
*
*  Notes:
*
*  1)  The proper motions in RA are dRA/dt rather than
*      cos(Dec)*dRA/dt, and are per year rather than per century.
*
*  2)  Note that conversion from Julian epoch 2000.0 to Besselian
*      epoch 1950.0 only is provided for.  Conversions involving
*      other epochs will require use of the appropriate precession,
*      proper motion, and E-terms routines before and/or after
*      FK524 is called.
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

      DOUBLE PRECISION R2000,D2000,DR2000,DD2000,P2000,V2000,
     :                 R1950,D1950,DR1950,DD1950,P1950,V1950


*  Miscellaneous
      DOUBLE PRECISION R,D,UR,UD,PX,RV
      DOUBLE PRECISION SR,CR,SD,CD,X,Y,Z,W
      DOUBLE PRECISION V1(6),V2(6)
      DOUBLE PRECISION XD,YD,ZD
      DOUBLE PRECISION RXYZ,WD,RXYSQ,RXY
      INTEGER I,J

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
      DOUBLE PRECISION A(6),EMI(6,6)
      DATA A/ -1.62557D-6,  -0.31919D-6, -0.13843D-6,
     :        +1.245D-3,    -1.580D-3,   -0.659D-3/

      DATA (EMI(I,1),I=1,6) / +0.9999256795D0,
     :                        -0.0111814828D0,
     :                        -0.0048590040D0,
     :                        -0.000551D0,
     :                        -0.238560D0,
     :                        +0.435730D0 /

      DATA (EMI(I,2),I=1,6) / +0.0111814828D0,
     :                        +0.9999374849D0,
     :                        -0.0000271557D0,
     :                        +0.238509D0,
     :                        -0.002667D0,
     :                        -0.008541D0 /

      DATA (EMI(I,3),I=1,6) / +0.0048590039D0,
     :                        -0.0000271771D0,
     :                        +0.9999881946D0,
     :                        -0.435614D0,
     :                        +0.012254D0,
     :                        +0.002117D0 /

      DATA (EMI(I,4),I=1,6) / -0.00000242389840D0,
     :                        +0.00000002710544D0,
     :                        +0.00000001177742D0,
     :                        +0.99990432D0,
     :                        -0.01118145D0,
     :                        -0.00485852D0 /

      DATA (EMI(I,5),I=1,6) / -0.00000002710544D0,
     :                        -0.00000242392702D0,
     :                        +0.00000000006585D0,
     :                        +0.01118145D0,
     :                        +0.99991613D0,
     :                        -0.00002716D0 /

      DATA (EMI(I,6),I=1,6) / -0.00000001177742D0,
     :                        +0.00000000006585D0,
     :                        -0.00000242404995D0,
     :                        +0.00485852D0,
     :                        -0.00002717D0,
     :                        +0.99996684D0 /



*  Pick up J2000 data (units radians and arcsec/JC)
      R=R2000
      D=D2000
      UR=DR2000*PMF
      UD=DD2000*PMF
      PX=P2000
      RV=V2000

*  Spherical to Cartesian
      SR=SIN(R)
      CR=COS(R)
      SD=SIN(D)
      CD=COS(D)

      X=CR*CD
      Y=SR*CD
      Z=   SD

      W=VF*RV*PX

      V1(1)=X
      V1(2)=Y
      V1(3)=Z

      V1(4)=-UR*Y-CR*SD*UD+W*X
      V1(5)= UR*X-SR*SD*UD+W*Y
      V1(6)=         CD*UD+W*Z

*  Convert position+velocity vector to BN system
      DO I=1,6
         W=0D0
         DO J=1,6
            W=W+EMI(I,J)*V1(J)
         END DO
         V2(I)=W
      END DO

*  Position vector components and magnitude
      X=V2(1)
      Y=V2(2)
      Z=V2(3)
      RXYZ=SQRT(X*X+Y*Y+Z*Z)

*  Apply E-terms to position
      W=X*A(1)+Y*A(2)+Z*A(3)
      X=X+A(1)*RXYZ-W*X
      Y=Y+A(2)*RXYZ-W*Y
      Z=Z+A(3)*RXYZ-W*Z

*  Recompute magnitude
      RXYZ=SQRT(X*X+Y*Y+Z*Z)

*  Apply E-terms to both position and velocity
      X=V2(1)
      Y=V2(2)
      Z=V2(3)
      W=X*A(1)+Y*A(2)+Z*A(3)
      WD=X*A(4)+Y*A(5)+Z*A(6)
      X=X+A(1)*RXYZ-W*X
      Y=Y+A(2)*RXYZ-W*Y
      Z=Z+A(3)*RXYZ-W*Z
      XD=V2(4)+A(4)*RXYZ-WD*X
      YD=V2(5)+A(5)*RXYZ-WD*Y
      ZD=V2(6)+A(6)*RXYZ-WD*Z

*  Convert to spherical
      RXYSQ=X*X+Y*Y
      RXY=SQRT(RXYSQ)

      IF (X.EQ.0D0.AND.Y.EQ.0D0) THEN
         R=0D0
      ELSE
         R=ATAN2(Y,X)
         IF (R.LT.0.0D0) R=R+D2PI
      END IF
      D=ATAN2(Z,RXY)

      IF (RXY.GT.TINY) THEN
         UR=(X*YD-Y*XD)/RXYSQ
         UD=(ZD*RXYSQ-Z*(X*XD+Y*YD))/((RXYSQ+Z*Z)*RXY)
      END IF

*  Radial velocity and parallax
      IF (PX.GT.TINY) THEN
         RV=(X*XD+Y*YD+Z*ZD)/(PX*VF*RXYZ)
         PX=PX/RXYZ
      END IF

*  Return results
      R1950=R
      D1950=D
      DR1950=UR/PMF
      DD1950=UD/PMF
      P1950=PX
      V1950=RV

      END
