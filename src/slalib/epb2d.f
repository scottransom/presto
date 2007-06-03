      DOUBLE PRECISION FUNCTION sla_EPB2D (EPB)
*+
*     - - - - - -
*      E P B 2 D
*     - - - - - -
*
*  Conversion of Besselian Epoch to Modified Julian Date
*  (double precision)
*
*  Given:
*     EPB      dp       Besselian Epoch
*
*  The result is the Modified Julian Date (JD - 2400000.5).
*
*  Reference:
*     Lieske,J.H., 1979. Astron.Astrophys.,73,282.
*
*  P.T.Wallace   Starlink   February 1984
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

      DOUBLE PRECISION EPB


      sla_EPB2D = 15019.81352D0 + (EPB-1900D0)*365.242198781D0

      END
