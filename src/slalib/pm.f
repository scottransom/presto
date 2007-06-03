      SUBROUTINE sla_PM (R0, D0, PR, PD, PX, RV, EP0, EP1, R1, D1)
*+
*     - - -
*      P M
*     - - -
*
*  Apply corrections for proper motion to a star RA,Dec
*  (double precision)
*
*  References:
*     1984 Astronomical Almanac, pp B39-B41.
*     (also Lederle & Schwan, Astron. Astrophys. 134,
*      1-6, 1984)
*
*  Given:
*     R0,D0    dp     RA,Dec at epoch EP0 (rad)
*     PR,PD    dp     proper motions:  RA,Dec changes per year of epoch
*     PX       dp     parallax (arcsec)
*     RV       dp     radial velocity (km/sec, +ve if receding)
*     EP0      dp     start epoch in years (e.g. Julian epoch)
*     EP1      dp     end epoch in years (same system as EP0)
*
*  Returned:
*     R1,D1    dp     RA,Dec at epoch EP1 (rad)
*
*  Called:
*     sla_DCS2C       spherical to Cartesian
*     sla_DCC2S       Cartesian to spherical
*     sla_DRANRM      normalize angle 0-2Pi
*
*  Notes:
*
*  1  The proper motions in RA are dRA/dt rather than cos(Dec)*dRA/dt,
*     and are in the same coordinate system as R0,D0.
*
*  2  If the available proper motions are pre-FK5 they will be per
*     tropical year rather than per Julian year, and so the epochs
*     must both be Besselian rather than Julian.  In such cases, a
*     scaling factor of 365.2422D0/365.25D0 should be applied to the
*     radial velocity before use.
*
*  P.T.Wallace   Starlink   19 January 2000
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

      DOUBLE PRECISION R0,D0,PR,PD,PX,RV,EP0,EP1,R1,D1

*  Km/s to AU/year multiplied by arcseconds to radians
      DOUBLE PRECISION VFR
      PARAMETER (VFR=(365.25D0*86400D0/149597870D0)*4.8481368111D-6)

      INTEGER I
      DOUBLE PRECISION sla_DRANRM
      DOUBLE PRECISION W,EM(3),T,P(3)



*  Spherical to Cartesian
      CALL sla_DCS2C(R0,D0,P)

*  Space motion (radians per year)
      W=VFR*RV*PX
      EM(1)=-PR*P(2)-PD*COS(R0)*SIN(D0)+W*P(1)
      EM(2)= PR*P(1)-PD*SIN(R0)*SIN(D0)+W*P(2)
      EM(3)=         PD*COS(D0)        +W*P(3)

*  Apply the motion
      T=EP1-EP0
      DO I=1,3
         P(I)=P(I)+T*EM(I)
      END DO

*  Cartesian to spherical
      CALL sla_DCC2S(P,R1,D1)
      R1=sla_DRANRM(R1)

      END
