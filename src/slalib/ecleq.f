      SUBROUTINE sla_ECLEQ (DL, DB, DATE, DR, DD)
*+
*     - - - - - -
*      E C L E Q
*     - - - - - -
*
*  Transformation from ecliptic coordinates to
*  J2000.0 equatorial coordinates (double precision)
*
*  Given:
*     DL,DB       dp      ecliptic longitude and latitude
*                           (mean of date, IAU 1980 theory, radians)
*     DATE        dp      TDB (loosely ET) as Modified Julian Date
*                                              (JD-2400000.5)
*  Returned:
*     DR,DD       dp      J2000.0 mean RA,Dec (radians)
*
*  Called:
*     sla_DCS2C, sla_ECMAT, sla_DIMXV, sla_PREC, sla_EPJ, sla_DCC2S,
*     sla_DRANRM, sla_DRANGE
*
*  P.T.Wallace   Starlink   March 1986
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

      DOUBLE PRECISION DL,DB,DATE,DR,DD

      DOUBLE PRECISION sla_EPJ,sla_DRANRM,sla_DRANGE

      DOUBLE PRECISION RMAT(3,3),V1(3),V2(3)



*  Spherical to Cartesian
      CALL sla_DCS2C(DL,DB,V1)

*  Ecliptic to equatorial
      CALL sla_ECMAT(DATE,RMAT)
      CALL sla_DIMXV(RMAT,V1,V2)

*  Mean of date to J2000
      CALL sla_PREC(2000D0,sla_EPJ(DATE),RMAT)
      CALL sla_DIMXV(RMAT,V2,V1)

*  Cartesian to spherical
      CALL sla_DCC2S(V1,DR,DD)

*  Express in conventional ranges
      DR=sla_DRANRM(DR)
      DD=sla_DRANGE(DD)

      END
