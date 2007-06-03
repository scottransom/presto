      SUBROUTINE sla_EQGAL (DR, DD, DL, DB)
*+
*     - - - - - -
*      E Q G A L
*     - - - - - -
*
*  Transformation from J2000.0 equatorial coordinates to
*  IAU 1958 galactic coordinates (double precision)
*
*  Given:
*     DR,DD       dp       J2000.0 RA,Dec
*
*  Returned:
*     DL,DB       dp       galactic longitude and latitude L2,B2
*
*  (all arguments are radians)
*
*  Called:
*     sla_DCS2C, sla_DMXV, sla_DCC2S, sla_DRANRM, sla_DRANGE
*
*  Note:
*     The equatorial coordinates are J2000.0.  Use the routine
*     sla_EG50 if conversion from B1950.0 'FK4' coordinates is
*     required.
*
*  Reference:
*     Blaauw et al, Mon.Not.R.Astron.Soc.,121,123 (1960)
*
*  P.T.Wallace   Starlink   21 September 1998
*
*  Copyright (C) 1998 Rutherford Appleton Laboratory
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

      DOUBLE PRECISION DR,DD,DL,DB

      DOUBLE PRECISION sla_DRANRM,sla_DRANGE

      DOUBLE PRECISION V1(3),V2(3)

*
*  L2,B2 system of galactic coordinates
*
*  P = 192.25       RA of galactic north pole (mean B1950.0)
*  Q =  62.6        inclination of galactic to mean B1950.0 equator
*  R =  33          longitude of ascending node
*
*  P,Q,R are degrees
*
*  Equatorial to galactic rotation matrix (J2000.0), obtained by
*  applying the standard FK4 to FK5 transformation, for zero proper
*  motion in FK5, to the columns of the B1950 equatorial to
*  galactic rotation matrix:
*
      DOUBLE PRECISION RMAT(3,3)
      DATA RMAT(1,1),RMAT(1,2),RMAT(1,3),
     :     RMAT(2,1),RMAT(2,2),RMAT(2,3),
     :     RMAT(3,1),RMAT(3,2),RMAT(3,3)/
     : -0.054875539726D0,-0.873437108010D0,-0.483834985808D0,
     : +0.494109453312D0,-0.444829589425D0,+0.746982251810D0,
     : -0.867666135858D0,-0.198076386122D0,+0.455983795705D0/



*  Spherical to Cartesian
      CALL sla_DCS2C(DR,DD,V1)

*  Equatorial to galactic
      CALL sla_DMXV(RMAT,V1,V2)

*  Cartesian to spherical
      CALL sla_DCC2S(V2,DL,DB)

*  Express in conventional ranges
      DL=sla_DRANRM(DL)
      DB=sla_DRANGE(DB)

      END
