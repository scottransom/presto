      SUBROUTINE sla__ATMT (R0, T0, ALPHA, GAMM2, DELM2,
     :                      C1, C2, C3, C4, C5, C6, R, T, DN, RDNDR)
*+
*     - - - - -
*      A T M T
*     - - - - -
*
*  Internal routine used by REFRO
*
*  Refractive index and derivative with respect to height for the
*  troposphere.
*
*  Given:
*    R0      d    height of observer from centre of the Earth (metre)
*    T0      d    temperature at the observer (K)
*    ALPHA   d    alpha          )
*    GAMM2   d    gamma minus 2  ) see HMNAO paper
*    DELM2   d    delta minus 2  )
*    C1      d    useful term  )
*    C2      d    useful term  )
*    C3      d    useful term  ) see source
*    C4      d    useful term  ) of sla_REFRO
*    C5      d    useful term  )
*    C6      d    useful term  )
*    R       d    current distance from the centre of the Earth (metre)
*
*  Returned:
*    T       d    temperature at R (K)
*    DN      d    refractive index at R
*    RDNDR   d    R * rate the refractive index is changing at R
*
*  Note that in the optical case C5 and C6 are zero.
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

      DOUBLE PRECISION R0,T0,ALPHA,GAMM2,DELM2,C1,C2,C3,C4,C5,C6,
     :                 R,T,DN,RDNDR

      DOUBLE PRECISION TT0,TT0GM2,TT0DM2


      T = MAX(MIN(T0-ALPHA*(R-R0),320D0),100D0)
      TT0 = T/T0
      TT0GM2 = TT0**GAMM2
      TT0DM2 = TT0**DELM2
      DN = 1D0+(C1*TT0GM2-(C2-C5/T)*TT0DM2)*TT0
      RDNDR = R*(-C3*TT0GM2+(C4-C6/TT0)*TT0DM2)

      END
