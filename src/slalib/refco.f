      SUBROUTINE sla_REFCO ( HM, TDK, PMB, RH, WL, PHI, TLR, EPS,
     :                       REFA, REFB )
*+
*     - - - - - -
*      R E F C O
*     - - - - - -
*
*  Determine the constants A and B in the atmospheric refraction
*  model dZ = A tan Z + B tan**3 Z.
*
*  Z is the "observed" zenith distance (i.e. affected by refraction)
*  and dZ is what to add to Z to give the "topocentric" (i.e. in vacuo)
*  zenith distance.
*
*  Given:
*    HM      d     height of the observer above sea level (metre)
*    TDK     d     ambient temperature at the observer (K)
*    PMB     d     pressure at the observer (millibar)
*    RH      d     relative humidity at the observer (range 0-1)
*    WL      d     effective wavelength of the source (micrometre)
*    PHI     d     latitude of the observer (radian, astronomical)
*    TLR     d     temperature lapse rate in the troposphere (K/metre)
*    EPS     d     precision required to terminate iteration (radian)
*
*  Returned:
*    REFA    d     tan Z coefficient (radian)
*    REFB    d     tan**3 Z coefficient (radian)
*
*  Called:  sla_REFRO
*
*  Notes:
*
*  1  Typical values for the TLR and EPS arguments might be 0.0065D0 and
*     1D-10 respectively.
*
*  2  The radio refraction is chosen by specifying WL > 100 micrometres.
*
*  3  The routine is a slower but more accurate alternative to the
*     sla_REFCOQ routine.  The constants it produces give perfect
*     agreement with sla_REFRO at zenith distances arctan(1) (45 deg)
*     and arctan(4) (about 76 deg).  It achieves 0.5 arcsec accuracy
*     for ZD < 80 deg, 0.01 arcsec accuracy for ZD < 60 deg, and
*     0.001 arcsec accuracy for ZD < 45 deg.
*
*  P.T.Wallace   Starlink   22 May 2004
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

      DOUBLE PRECISION HM,TDK,PMB,RH,WL,PHI,TLR,EPS,REFA,REFB

      DOUBLE PRECISION ATN1,ATN4,R1,R2

*  Sample zenith distances: arctan(1) and arctan(4)
      PARAMETER (ATN1=0.7853981633974483D0,
     :           ATN4=1.325817663668033D0)



*  Determine refraction for the two sample zenith distances
      CALL sla_REFRO(ATN1,HM,TDK,PMB,RH,WL,PHI,TLR,EPS,R1)
      CALL sla_REFRO(ATN4,HM,TDK,PMB,RH,WL,PHI,TLR,EPS,R2)

*  Solve for refraction constants
      REFA = (64D0*R1-R2)/60D0
      REFB = (R2-4D0*R1)/60D0

      END
