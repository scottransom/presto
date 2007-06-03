      DOUBLE PRECISION FUNCTION sla_DT (EPOCH)
*+
*     - - -
*      D T
*     - - -
*
*  Estimate the offset between dynamical time and Universal Time
*  for a given historical epoch.
*
*  Given:
*     EPOCH       d        (Julian) epoch (e.g. 1850D0)
*
*  The result is a rough estimate of ET-UT (after 1984, TT-UT) at
*  the given epoch, in seconds.
*
*  Notes:
*
*  1  Depending on the epoch, one of three parabolic approximations
*     is used:
*
*      before 979    Stephenson & Morrison's 390 BC to AD 948 model
*      979 to 1708   Stephenson & Morrison's 948 to 1600 model
*      after 1708    McCarthy & Babcock's post-1650 model
*
*     The breakpoints are chosen to ensure continuity:  they occur
*     at places where the adjacent models give the same answer as
*     each other.
*
*  2  The accuracy is modest, with errors of up to 20 sec during
*     the interval since 1650, rising to perhaps 30 min by 1000 BC.
*     Comparatively accurate values from AD 1600 are tabulated in
*     the Astronomical Almanac (see section K8 of the 1995 AA).
*
*  3  The use of double-precision for both argument and result is
*     purely for compatibility with other SLALIB time routines.
*
*  4  The models used are based on a lunar tidal acceleration value
*     of -26.00 arcsec per century.
*
*  Reference:  Explanatory Supplement to the Astronomical Almanac,
*              ed P.K.Seidelmann, University Science Books (1992),
*              section 2.553, p83.  This contains references to
*              the Stephenson & Morrison and McCarthy & Babcock
*              papers.
*
*  P.T.Wallace   Starlink   1 March 1995
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

      DOUBLE PRECISION EPOCH
      DOUBLE PRECISION T,W,S


*  Centuries since 1800
      T=(EPOCH-1800D0)/100D0

*  Select model
      IF (EPOCH.GE.1708.185161980887D0) THEN

*     Post-1708: use McCarthy & Babcock
         W=T-0.19D0
         S=5.156D0+13.3066D0*W*W
      ELSE IF (EPOCH.GE.979.0258204760233D0) THEN

*     979-1708: use Stephenson & Morrison's 948-1600 model
         S=25.5D0*T*T
      ELSE

*     Pre-979: use Stephenson & Morrison's 390 BC to AD 948 model
         S=1360.0D0+(320D0+44.3D0*T)*T
      END IF

*  Result
      sla_DT=S

      END
