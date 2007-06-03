      SUBROUTINE sla_GALSUP (DL, DB, DSL, DSB)
*+
*     - - - - - - -
*      G A L S U P
*     - - - - - - -
*
*  Transformation from IAU 1958 galactic coordinates to
*  de Vaucouleurs supergalactic coordinates (double precision)
*
*  Given:
*     DL,DB       dp       galactic longitude and latitude L2,B2
*
*  Returned:
*     DSL,DSB     dp       supergalactic longitude and latitude
*
*  (all arguments are radians)
*
*  Called:
*     sla_DCS2C, sla_DMXV, sla_DCC2S, sla_DRANRM, sla_DRANGE
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
*  P.T.Wallace   Starlink   25 January 1999
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

      DOUBLE PRECISION DL,DB,DSL,DSB

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
      CALL sla_DCS2C(DL,DB,V1)

*  Galactic to supergalactic
      CALL sla_DMXV(RMAT,V1,V2)

*  Cartesian to spherical
      CALL sla_DCC2S(V2,DSL,DSB)

*  Express in conventional ranges
      DSL=sla_DRANRM(DSL)
      DSB=sla_DRANGE(DSB)

      END
