      SUBROUTINE sla__ATMS (RT, TT, DNT, GAMAL, R, DN, RDNDR)
*+
*     - - - - -
*      A T M S
*     - - - - -
*
*  Internal routine used by REFRO
*
*  Refractive index and derivative with respect to height for the
*  stratosphere.
*
*  Given:
*    RT      d    height of tropopause from centre of the Earth (metre)
*    TT      d    temperature at the tropopause (K)
*    DNT     d    refractive index at the tropopause
*    GAMAL   d    constant of the atmospheric model = G*MD/R
*    R       d    current distance from the centre of the Earth (metre)
*
*  Returned:
*    DN      d    refractive index at R
*    RDNDR   d    R * rate the refractive index is changing at R
*
*  Last revision:   26 December 2004
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

      DOUBLE PRECISION RT,TT,DNT,GAMAL,R,DN,RDNDR

      DOUBLE PRECISION B,W


      B = GAMAL/TT
      W = (DNT-1D0)*EXP(-B*(R-RT))
      DN = 1D0+W
      RDNDR = -R*B*W

      END
