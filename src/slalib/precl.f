      SUBROUTINE sla_PRECL (EP0, EP1, RMATP)
*+
*     - - - - - -
*      P R E C L
*     - - - - - -
*
*  Form the matrix of precession between two epochs, using the
*  model of Simon et al (1994), which is suitable for long
*  periods of time.
*
*  (double precision)
*
*  Given:
*     EP0    dp         beginning epoch
*     EP1    dp         ending epoch
*
*  Returned:
*     RMATP  dp(3,3)    precession matrix
*
*  Notes:
*
*     1)  The epochs are TDB Julian epochs.
*
*     2)  The matrix is in the sense   V(EP1)  =  RMATP * V(EP0)
*
*     3)  The absolute accuracy of the model is limited by the
*         uncertainty in the general precession, about 0.3 arcsec per
*         1000 years.  The remainder of the formulation provides a
*         precision of 1 mas over the interval from 1000AD to 3000AD,
*         0.1 arcsec from 1000BC to 5000AD and 1 arcsec from
*         4000BC to 8000AD.
*
*  Reference:
*     Simon, J.L. et al., 1994. Astron.Astrophys., 282, 663-683.
*
*  Called:  sla_DEULER
*
*  P.T.Wallace   Starlink   23 August 1996
*
*  Copyright (C) 1996 Rutherford Appleton Laboratory
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

      DOUBLE PRECISION EP0,EP1,RMATP(3,3)

*  Arc seconds to radians
      DOUBLE PRECISION AS2R
      PARAMETER (AS2R=0.484813681109535994D-5)

      DOUBLE PRECISION T0,T,TAS2R,W,ZETA,Z,THETA



*  Interval between basic epoch J2000.0 and beginning epoch (1000JY)
      T0 = (EP0-2000D0)/1000D0

*  Interval over which precession required (1000JY)
      T = (EP1-EP0)/1000D0

*  Euler angles
      TAS2R = T*AS2R
      W =      23060.9097D0+
     :          (139.7459D0+
     :           (-0.0038D0+
     :           (-0.5918D0+
     :           (-0.0037D0+
     :             0.0007D0*T0)*T0)*T0)*T0)*T0

      ZETA =   (W+(30.2226D0+
     :            (-0.2523D0+
     :            (-0.3840D0+
     :            (-0.0014D0+
     :              0.0007D0*T0)*T0)*T0)*T0+
     :            (18.0183D0+
     :            (-0.1326D0+
     :             (0.0006D0+
     :              0.0005D0*T0)*T0)*T0+
     :            (-0.0583D0+
     :            (-0.0001D0+
     :              0.0007D0*T0)*T0+
     :            (-0.0285D0+
     :            (-0.0002D0)*T)*T)*T)*T)*T)*TAS2R

      Z =     (W+(109.5270D0+
     :             (0.2446D0+
     :            (-1.3913D0+
     :            (-0.0134D0+
     :              0.0026D0*T0)*T0)*T0)*T0+
     :            (18.2667D0+
     :            (-1.1400D0+
     :            (-0.0173D0+
     :              0.0044D0*T0)*T0)*T0+
     :            (-0.2821D0+
     :            (-0.0093D0+
     :              0.0032D0*T0)*T0+
     :            (-0.0301D0+
     :              0.0006D0*T0
     :             -0.0001D0*T)*T)*T)*T)*T)*TAS2R

      THETA =  (20042.0207D0+
     :           (-85.3131D0+
     :            (-0.2111D0+
     :             (0.3642D0+
     :             (0.0008D0+
     :            (-0.0005D0)*T0)*T0)*T0)*T0)*T0+
     :           (-42.6566D0+
     :            (-0.2111D0+
     :             (0.5463D0+
     :             (0.0017D0+
     :            (-0.0012D0)*T0)*T0)*T0)*T0+
     :           (-41.8238D0+
     :             (0.0359D0+
     :             (0.0027D0+
     :            (-0.0001D0)*T0)*T0)*T0+
     :            (-0.0731D0+
     :             (0.0019D0+
     :              0.0009D0*T0)*T0+
     :            (-0.0127D0+
     :              0.0011D0*T0+0.0004D0*T)*T)*T)*T)*T)*TAS2R

*  Rotation matrix
      CALL sla_DEULER('ZYZ',-ZETA,THETA,-Z,RMATP)

      END
