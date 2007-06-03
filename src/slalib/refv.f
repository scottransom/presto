      SUBROUTINE sla_REFV (VU, REFA, REFB, VR)
*+
*     - - - - -
*      R E F V
*     - - - - -
*
*  Adjust an unrefracted Cartesian vector to include the effect of
*  atmospheric refraction, using the simple A tan Z + B tan**3 Z
*  model.
*
*  Given:
*    VU    dp    unrefracted position of the source (Az/El 3-vector)
*    REFA  dp    tan Z coefficient (radian)
*    REFB  dp    tan**3 Z coefficient (radian)
*
*  Returned:
*    VR    dp    refracted position of the source (Az/El 3-vector)
*
*  Notes:
*
*  1  This routine applies the adjustment for refraction in the
*     opposite sense to the usual one - it takes an unrefracted
*     (in vacuo) position and produces an observed (refracted)
*     position, whereas the A tan Z + B tan**3 Z model strictly
*     applies to the case where an observed position is to have the
*     refraction removed.  The unrefracted to refracted case is
*     harder, and requires an inverted form of the text-book
*     refraction models;  the algorithm used here is equivalent to
*     one iteration of the Newton-Raphson method applied to the above
*     formula.
*
*  2  Though optimized for speed rather than precision, the present
*     routine achieves consistency with the refracted-to-unrefracted
*     A tan Z + B tan**3 Z model at better than 1 microarcsecond within
*     30 degrees of the zenith and remains within 1 milliarcsecond to
*     beyond ZD 70 degrees.  The inherent accuracy of the model is, of
*     course, far worse than this - see the documentation for sla_REFCO
*     for more information.
*
*  3  At low elevations (below about 3 degrees) the refraction
*     correction is held back to prevent arithmetic problems and
*     wildly wrong results.  For optical/IR wavelengths, over a wide
*     range of observer heights and corresponding temperatures and
*     pressures, the following levels of accuracy (arcsec, worst case)
*     are achieved, relative to numerical integration through a model
*     atmosphere:
*
*              ZD    error
*
*              80      0.7
*              81      1.3
*              82      2.5
*              83      5
*              84     10
*              85     20
*              86     55
*              87    160
*              88    360
*              89    640
*              90   1100
*              91   1700         } relevant only to
*              92   2600         } high-elevation sites
*
*     The results for radio are slightly worse over most of the range,
*     becoming significantly worse below ZD=88 and unusable beyond
*     ZD=90.
*
*  4  See also the routine sla_REFZ, which performs the adjustment to
*     the zenith distance rather than in Cartesian Az/El coordinates.
*     The present routine is faster than sla_REFZ and, except very low down,
*     is equally accurate for all practical purposes.  However, beyond
*     about ZD 84 degrees sla_REFZ should be used, and for the utmost
*     accuracy iterative use of sla_REFRO should be considered.
*
*  P.T.Wallace   Starlink   10 April 2004
*
*  Copyright (C) 2004 Rutherford Appleton Laboratory
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

      DOUBLE PRECISION VU(3),REFA,REFB,VR(3)

      DOUBLE PRECISION X,Y,Z1,Z,ZSQ,RSQ,R,WB,WT,D,CD,F



*  Initial estimate = unrefracted vector
      X = VU(1)
      Y = VU(2)
      Z1 = VU(3)

*  Keep correction approximately constant below about 3 deg elevation
      Z = MAX(Z1,0.05D0)

*  One Newton-Raphson iteration
      ZSQ = Z*Z
      RSQ = X*X+Y*Y
      R = SQRT(RSQ)
      WB = REFB*RSQ/ZSQ
      WT = (REFA+WB)/(1D0+(REFA+3D0*WB)*(ZSQ+RSQ)/ZSQ)
      D = WT*R/Z
      CD = 1D0-D*D/2D0
      F = CD*(1D0-WT)

*  Post-refraction x,y,z
      VR(1) = X*F
      VR(2) = Y*F
      VR(3) = CD*(Z+D*R)+(Z1-Z)

      END
