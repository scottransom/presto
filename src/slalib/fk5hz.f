      SUBROUTINE sla_FK5HZ (R5,D5,EPOCH,RH,DH)
*+
*     - - - - - -
*      F K 5 H Z
*     - - - - - -
*
*  Transform an FK5 (J2000) star position into the frame of the
*  Hipparcos catalogue, assuming zero Hipparcos proper motion.
*
*  (double precision)
*
*  This routine converts a star position from the FK5 system to
*  the Hipparcos system, in such a way that the Hipparcos proper
*  motion is zero.  Because such a star has, in general, a non-zero
*  proper motion in the FK5 system, the routine requires the epoch
*  at which the position in the FK5 system was determined.
*
*  Given:
*     R5        d      FK5 RA (radians), equinox J2000, epoch EPOCH
*     D5        d      FK5 Dec (radians), equinox J2000, epoch EPOCH
*     EPOCH     d      Julian epoch (TDB)
*
*  Returned (all Hipparcos):
*     RH        d      RA (radians)
*     DH        d      Dec (radians)
*
*  Called:  sla_DCS2C, sla_DAV2M, sla_DIMXV, sla_DMXV, sla_DCC2S,
*           sla_DRANRM
*
*  Notes:
*
*  1)  The FK5 to Hipparcos transformation consists of a pure
*      rotation and spin;  zonal errors in the FK5 catalogue are
*      not taken into account.
*
*  2)  The published orientation and spin components are interpreted
*      as "axial vectors".  An axial vector points at the pole of the
*      rotation and its length is the amount of rotation in radians.
*
*  3)  See also sla_FK52H, sla_H2FK5, sla_HFK5Z.
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

      DOUBLE PRECISION R5,D5,EPOCH,RH,DH

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

      DOUBLE PRECISION P5E(3),ORTN(3),R5H(3,3),T,VST(3),RST(3,3),
     :                 P5(3),PH(3),W

      DOUBLE PRECISION sla_DRANRM



*  FK5 barycentric position vector.
      CALL sla_DCS2C(R5,D5,P5E)

*  FK5 to Hipparcos orientation matrix.
      ORTN(1) = EPX
      ORTN(2) = EPY
      ORTN(3) = EPZ
      CALL sla_DAV2M(ORTN,R5H)

*  Time interval from epoch to J2000.
      T = 2000D0-EPOCH

*  Axial vector:  accumulated Hipparcos wrt FK5 spin over that interval.
      VST(1) = OMX*T
      VST(2) = OMY*T
      VST(3) = OMZ*T

*  Express the accumulated spin as a rotation matrix.
      CALL sla_DAV2M(VST,RST)

*  Derotate the vector's FK5 axes back to epoch.
      CALL sla_DIMXV(RST,P5E,P5)

*  Rotate the vector into the Hipparcos frame.
      CALL sla_DMXV(R5H,P5,PH)

*  Hipparcos vector to spherical.
      CALL sla_DCC2S(PH,W,DH)
      RH = sla_DRANRM(W)

      END
