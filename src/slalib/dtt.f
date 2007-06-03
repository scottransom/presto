      DOUBLE PRECISION FUNCTION sla_DTT (UTC)
*+
*     - - - -
*      D T T
*     - - - -
*
*  Increment to be applied to Coordinated Universal Time UTC to give
*  Terrestrial Time TT (formerly Ephemeris Time ET)
*
*  (double precision)
*
*  Given:
*     UTC      d      UTC date as a modified JD (JD-2400000.5)
*
*  Result:  TT-UTC in seconds
*
*  Notes:
*
*  1  The UTC is specified to be a date rather than a time to indicate
*     that care needs to be taken not to specify an instant which lies
*     within a leap second.  Though in most cases UTC can include the
*     fractional part, correct behaviour on the day of a leap second
*     can only be guaranteed up to the end of the second 23:59:59.
*
*  2  Pre 1972 January 1 a fixed value of 10 + ET-TAI is returned.
*
*  3  See also the routine sla_DT, which roughly estimates ET-UT for
*     historical epochs.
*
*  Called:  sla_DAT
*
*  P.T.Wallace   Starlink   6 December 1994
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

      DOUBLE PRECISION UTC

      DOUBLE PRECISION sla_DAT


      sla_DTT=32.184D0+sla_DAT(UTC)

      END
