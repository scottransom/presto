      SUBROUTINE sla_SUPGAL (DSL, DSB, DL, DB)
*+
*     - - - - - - -
*      S U P G A L
*     - - - - - - -
*
*  Transformation from de Vaucouleurs supergalactic coordinates
*  to IAU 1958 galactic coordinates (double precision)
*
*  Given:
*     DSL,DSB     dp       supergalactic longitude and latitude
*
*  Returned:
*     DL,DB       dp       galactic longitude and latitude L2,B2
*
*  (all arguments are radians)
*
*  Called:
*     sla_DCS2C, sla_DIMXV, sla_DCC2S, sla_DRANRM, sla_DRANGE
*
*  References:
*
*     de Vaucouleurs, de Vaucouleurs, & Corwin, Second Reference
*     Catalogue of Bright Galaxies, U. Texas, page 8.
*
*     Systems & Applied Sciences Corp., Documentation for the
*     machine-readable version of the above catalogue,
*     Contract NAS 5-26490.
*
*    (These two references give different values for the galactic
*     longitude of the supergalactic origin.  Both are wrong;  the
*     correct value is L2=137.37.)
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

      DOUBLE PRECISION DSL,DSB,DL,DB

      DOUBLE PRECISION sla_DRANRM,sla_DRANGE

      DOUBLE PRECISION V1(3),V2(3)

*
*  System of supergalactic coordinates:
*
*    SGL   SGB        L2     B2      (deg)
*     -    +90      47.37  +6.32
*     0     0         -      0
*
*  Galactic to supergalactic rotation matrix:
*
      DOUBLE PRECISION RMAT(3,3)
      DATA RMAT(1,1),RMAT(1,2),RMAT(1,3),
     :     RMAT(2,1),RMAT(2,2),RMAT(2,3),
     :     RMAT(3,1),RMAT(3,2),RMAT(3,3)/
     : -0.735742574804D0,+0.677261296414D0,+0.000000000000D0,
     : -0.074553778365D0,-0.080991471307D0,+0.993922590400D0,
     : +0.673145302109D0,+0.731271165817D0,+0.110081262225D0/



*  Spherical to Cartesian
      CALL sla_DCS2C(DSL,DSB,V1)

*  Supergalactic to galactic
      CALL sla_DIMXV(RMAT,V1,V2)

*  Cartesian to spherical
      CALL sla_DCC2S(V2,DL,DB)

*  Express in conventional ranges
      DL=sla_DRANRM(DL)
      DB=sla_DRANGE(DB)

      END
