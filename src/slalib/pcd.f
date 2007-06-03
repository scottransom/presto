      SUBROUTINE sla_PCD (DISCO,X,Y)
*+
*     - - - -
*      P C D
*     - - - -
*
*  Apply pincushion/barrel distortion to a tangent-plane [x,y].
*
*  Given:
*     DISCO    d      pincushion/barrel distortion coefficient
*     X,Y      d      tangent-plane coordinates
*
*  Returned:
*     X,Y      d      distorted coordinates
*
*  Notes:
*
*  1)  The distortion is of the form RP = R*(1 + C*R**2), where R is
*      the radial distance from the tangent point, C is the DISCO
*      argument, and RP is the radial distance in the presence of
*      the distortion.
*
*  2)  For pincushion distortion, C is +ve;  for barrel distortion,
*      C is -ve.
*
*  3)  For X,Y in units of one projection radius (in the case of
*      a photographic plate, the focal length), the following
*      DISCO values apply:
*
*          Geometry          DISCO
*
*          astrograph         0.0
*          Schmidt           -0.3333
*          AAT PF doublet  +147.069
*          AAT PF triplet  +178.585
*          AAT f/8          +21.20
*          JKT f/8          +13.32
*
*  4)  There is a companion routine, sla_UNPCD, which performs the
*      inverse operation.
*
*  P.T.Wallace   Starlink   3 September 2000
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

      DOUBLE PRECISION DISCO,X,Y

      DOUBLE PRECISION F



      F=1D0+DISCO*(X*X+Y*Y)
      X=X*F
      Y=Y*F

      END
