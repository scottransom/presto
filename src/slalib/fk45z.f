      SUBROUTINE sla_FK45Z (R1950,D1950,BEPOCH,R2000,D2000)
*+
*     - - - - - -
*      F K 4 5 Z
*     - - - - - -
*
*  Convert B1950.0 FK4 star data to J2000.0 FK5 assuming zero
*  proper motion in the FK5 frame (double precision)
*
*  This routine converts stars from the old, Bessel-Newcomb, FK4
*  system to the new, IAU 1976, FK5, Fricke system, in such a
*  way that the FK5 proper motion is zero.  Because such a star
*  has, in general, a non-zero proper motion in the FK4 system,
*  the routine requires the epoch at which the position in the
*  FK4 system was determined.
*
*  The method is from Appendix 2 of Ref 1, but using the constants
*  of Ref 4.
*
*  Given:
*     R1950,D1950     dp    B1950.0 FK4 RA,Dec at epoch (rad)
*     BEPOCH          dp    Besselian epoch (e.g. 1979.3D0)
*
*  Returned:
*     R2000,D2000     dp    J2000.0 FK5 RA,Dec (rad)
*
*  Notes:
*
*  1)  The epoch BEPOCH is strictly speaking Besselian, but
*      if a Julian epoch is supplied the result will be
*      affected only to a negligible extent.
*
*  2)  Conversion from Besselian epoch 1950.0 to Julian epoch
*      2000.0 only is provided for.  Conversions involving other
*      epochs will require use of the appropriate precession,
*      proper motion, and E-terms routines before and/or
*      after FK45Z is called.
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
*     1  Aoki,S., et al, 1983.  Astron.Astrophys., 128, 263.
*
*     2  Smith, C.A. et al, 1989.  "The transformation of astrometric
*        catalog systems to the equinox J2000.0".  Astron.J. 97, 265.
*
*     3  Yallop, B.D. et al, 1989.  "Transformation of mean star places
*        from FK4 B1950.0 to FK5 J2000.0 using matrices in 6-space".
*        Astron.J. 97, 274.
*
*     4  Seidelmann, P.K. (ed), 1992.  "Explanatory Supplement to
*        the Astronomical Almanac", ISBN 0-935702-68-7.
*
*  Called:  sla_DCS2C, sla_EPJ, sla_EPB2D, sla_DCC2S, sla_DRANRM
*
*  P.T.Wallace   Starlink   21 September 1998
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

      DOUBLE PRECISION R1950,D1950,BEPOCH,R2000,D2000

      DOUBLE PRECISION D2PI
      PARAMETER (D2PI=6.283185307179586476925287D0)

      DOUBLE PRECISION W
      INTEGER I,J

*  Position and position+velocity vectors
      DOUBLE PRECISION R0(3),A1(3),V1(3),V2(6)

*  Radians per year to arcsec per century
      DOUBLE PRECISION PMF
      PARAMETER (PMF=100D0*60D0*60D0*360D0/D2PI)

*  Functions
      DOUBLE PRECISION sla_EPJ,sla_EPB2D,sla_DRANRM

*
*  CANONICAL CONSTANTS  (see references)
*

*  Vectors A and Adot, and matrix M (only half of which is needed here)
      DOUBLE PRECISION A(3),AD(3),EM(6,3)
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



*  Spherical to Cartesian
      CALL sla_DCS2C(R1950,D1950,R0)

*  Adjust vector A to give zero proper motion in FK5
      W=(BEPOCH-1950D0)/PMF
      DO I=1,3
         A1(I)=A(I)+W*AD(I)
      END DO

*  Remove e-terms
      W=R0(1)*A1(1)+R0(2)*A1(2)+R0(3)*A1(3)
      DO I=1,3
         V1(I)=R0(I)-A1(I)+W*R0(I)
      END DO

*  Convert position vector to Fricke system
      DO I=1,6
         W=0D0
         DO J=1,3
            W=W+EM(I,J)*V1(J)
         END DO
         V2(I)=W
      END DO

*  Allow for fictitious proper motion in FK4
      W=(sla_EPJ(sla_EPB2D(BEPOCH))-2000D0)/PMF
      DO I=1,3
         V2(I)=V2(I)+W*V2(I+3)
      END DO

*  Revert to spherical coordinates
      CALL sla_DCC2S(V2,W,D2000)
      R2000=sla_DRANRM(W)

      END
