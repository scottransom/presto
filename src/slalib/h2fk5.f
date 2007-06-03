      SUBROUTINE sla_H2FK5 (RH,DH,DRH,DDH,R5,D5,DR5,DD5)
*+
*     - - - - - -
*      H 2 F K 5
*     - - - - - -
*
*  Transform Hipparcos star data into the FK5 (J2000) system.
*
*  (double precision)
*
*  This routine transforms Hipparcos star positions and proper
*  motions into FK5 J2000.
*
*  Given (all Hipparcos, epoch J2000):
*     RH        d      RA (radians)
*     DH        d      Dec (radians)
*     DRH       d      proper motion in RA (dRA/dt, rad/Jyear)
*     DDH       d      proper motion in Dec (dDec/dt, rad/Jyear)
*
*  Returned (all FK5, equinox J2000, epoch J2000):
*     R5        d      RA (radians)
*     D5        d      Dec (radians)
*     DR5       d      proper motion in RA (dRA/dt, rad/Jyear)
*     DD5       d      proper motion in Dec (dDec/dt, rad/Jyear)
*
*  Called:  sla_DS2C6, sla_DAV2M, sla_DMXV, sla_DIMXV, sla_DVXV,
*           sla_DC62S, sla_DRANRM
*
*  Notes:
*
*  1)  The proper motions in RA are dRA/dt rather than
*      cos(Dec)*dRA/dt, and are per year rather than per century.
*
*  2)  The FK5 to Hipparcos transformation consists of a pure
*      rotation and spin;  zonal errors in the FK5 catalogue are
*      not taken into account.
*
*  3)  The published orientation and spin components are interpreted
*      as "axial vectors".  An axial vector points at the pole of the
*      rotation and its length is the amount of rotation in radians.
*
*  4)  See also sla_FK52H, sla_FK5HZ, sla_HFK5Z.
*
*  Reference:
*
*     M.Feissel & F.Mignard, Astron. Astrophys. 331, L33-L36 (1998).
*
*  P.T.Wallace   Starlink   22 June 1999
*
*  Copyright (C) 1999 Rutherford Appleton Laboratory
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

      DOUBLE PRECISION RH,DH,DRH,DDH,R5,D5,DR5,DD5

      DOUBLE PRECISION AS2R
      PARAMETER (AS2R=0.484813681109535994D-5)

*  FK5 to Hipparcos orientation and spin (radians, radians/year)
      DOUBLE PRECISION EPX,EPY,EPZ
      DOUBLE PRECISION OMX,OMY,OMZ

      PARAMETER ( EPX = -19.9D-3 * AS2R,
     :            EPY =  -9.1D-3 * AS2R,
     :            EPZ = +22.9D-3 * AS2R )

      PARAMETER ( OMX = -0.30D-3 * AS2R,
     :            OMY = +0.60D-3 * AS2R,
     :            OMZ = +0.70D-3 * AS2R )

      DOUBLE PRECISION PVH(6),ORTN(3),R5H(3,3),S5(3),SH(3),VV(3),
     :                 PV5(6),W,R,V
      INTEGER I

      DOUBLE PRECISION sla_DRANRM



*  Hipparcos barycentric position/velocity 6-vector (normalized).
      CALL sla_DS2C6(RH,DH,1D0,DRH,DDH,0D0,PVH)

*  FK5 to Hipparcos orientation matrix.
      ORTN(1) = EPX
      ORTN(2) = EPY
      ORTN(3) = EPZ
      CALL sla_DAV2M(ORTN,R5H)

*  Hipparcos wrt FK5 spin vector.
      S5(1) = OMX
      S5(2) = OMY
      S5(3) = OMZ

*  Rotate the spin vector into the Hipparcos frame.
      CALL sla_DMXV(R5H,S5,SH)

*  De-orient & de-spin the 6-vector into FK5 J2000.
      CALL sla_DIMXV(R5H,PVH,PV5)
      CALL sla_DVXV(PVH,SH,VV)
      DO I=1,3
         VV(I) = PVH(I+3)-VV(I)
      END DO
      CALL sla_DIMXV(R5H,VV,PV5(4))

*  FK5 6-vector to spherical.
      CALL sla_DC62S(PV5,W,D5,R,DR5,DD5,V)
      R5 = sla_DRANRM(W)

      END
