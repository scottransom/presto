      SUBROUTINE sla_WAIT (DELAY)
*+
*     - - - - -
*      W A I T
*     - - - - -
*
*  Interval wait
*
*  !!! Version for: SPARC/SunOS4, 
*                   SPARC/Solaris2, 
*                   DEC Mips/Ultrix
*                   DEC AXP/Digital Unix
*                   Intel/Linux
*                   Convex
*
*  Given:
*     DELAY     real      delay in seconds
*
*  Called:  SLEEP (a Fortran Intrinsic on all obove platforms)
*
*  P.T.Wallace   Starlink   22 January 1998
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

      REAL DELAY

      CALL SLEEP(NINT(DELAY))

      END
