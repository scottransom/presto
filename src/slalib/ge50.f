      SUBROUTINE sla_GE50 (DL, DB, DR, DD)
*+
*     - - - - -
*      G E 5 0
*     - - - - -
*
*  Transformation from IAU 1958 galactic coordinates to
*  B1950.0 'FK4' equatorial coordinates (double precision)
*
*  Given:
*     DL,DB       dp       galactic longitude and latitude L2,B2
*
*  Returned:
*     DR,DD       dp       B1950.0 'FK4' RA,Dec
*
*  (all arguments are radians)
*
*  Called:
*     sla_DCS2C, sla_DIMXV, sla_DCC2S, sla_ADDET, sla_DRANRM, sla_DRANGE
*
*  Note:
*     The equatorial coordinates are B1950.0 'FK4'.  Use the
*     routine sla_GALEQ if conversion to J2000.0 coordinates
*     is required.
*
*  Reference:
*     Blaauw et al, Mon.Not.R.Astron.Soc.,121,123 (1960)
*
*  P.T.Wallace   Starlink   5 September 1993
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

      DOUBLE PRECISION DL,DB,DR,DD

      DOUBLE PRECISION sla_DRANRM,sla_DRANGE

      DOUBLE PRECISION V1(3),V2(3),R,D,RE,DE

*
*  L2,B2 system of galactic coordinates
*
*  P = 192.25       RA of galactic north pole (mean B1950.0)
*  Q =  62.6        inclination of galactic to mean B1950.0 equator
*  R =  33          longitude of ascending node
*
*  P,Q,R are degrees
*
*
*  Equatorial to galactic rotation matrix
*
*  The Euler angles are P, Q, 90-R, about the z then y then
*  z axes.
*
*         +CP.CQ.SR-SP.CR     +SP.CQ.SR+CP.CR     -SQ.SR
*
*         -CP.CQ.CR-SP.SR     -SP.CQ.CR+CP.SR     +SQ.CR
*
*         +CP.SQ              +SP.SQ              +CQ
*

      DOUBLE PRECISION RMAT(3,3)
      DATA RMAT(1,1),RMAT(1,2),RMAT(1,3),
     :     RMAT(2,1),RMAT(2,2),RMAT(2,3),
     :     RMAT(3,1),RMAT(3,2),RMAT(3,3) /
     : -0.066988739415D0,-0.872755765852D0,-0.483538914632D0,
     : +0.492728466075D0,-0.450346958020D0,+0.744584633283D0,
     : -0.867600811151D0,-0.188374601723D0,+0.460199784784D0 /



*  Spherical to Cartesian
      CALL sla_DCS2C(DL,DB,V1)

*  Rotate to mean B1950.0
      CALL sla_DIMXV(RMAT,V1,V2)

*  Cartesian to spherical
      CALL sla_DCC2S(V2,R,D)

*  Introduce E-terms
      CALL sla_ADDET(R,D,1950D0,RE,DE)

*  Express in conventional ranges
      DR=sla_DRANRM(RE)
      DD=sla_DRANGE(DE)

      END
