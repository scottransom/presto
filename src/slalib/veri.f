      INTEGER FUNCTION sla_VERI ()
*+
*     - - - - -
*      V E R I
*     - - - - -
*
*  Report the SLALIB version number as an integer.
*
*  Given:
*    None
*
*  The result is the SLALIB version number as an integer m*1e6+n*1e3+r,
*  where m is the major version, n the minor version and r the release
*  number.
*
*  Notes:
*
*    To obtain the version number in a printable form, see
*    subroutine sla_vers(version).
*
*    The sla_veri subroutine was introduced in SLALIB version 2.5-1, so
*    if this function is absent, one can only tell that the release
*    predates that one.
*
*  Norman Gray   Starlink   8 April 2005
*
*  Copyright (C) 2005 Council for the Central Laboratory of the
*  Research Councils
*
*  Licence:
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

      sla_VERI=2005004

      END
