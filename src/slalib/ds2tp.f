      SUBROUTINE sla_DS2TP (RA, DEC, RAZ, DECZ, XI, ETA, J)
*+
*     - - - - - -
*      D S 2 T P
*     - - - - - -
*
*  Projection of spherical coordinates onto tangent plane:
*  "gnomonic" projection - "standard coordinates" (double precision)
*
*  Given:
*     RA,DEC      dp   spherical coordinates of point to be projected
*     RAZ,DECZ    dp   spherical coordinates of tangent point
*
*  Returned:
*     XI,ETA      dp   rectangular coordinates on tangent plane
*     J           int  status:   0 = OK, star on tangent plane
*                                1 = error, star too far from axis
*                                2 = error, antistar on tangent plane
*                                3 = error, antistar too far from axis
*
*  P.T.Wallace   Starlink   18 July 1996
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

      DOUBLE PRECISION RA,DEC,RAZ,DECZ,XI,ETA
      INTEGER J

      DOUBLE PRECISION SDECZ,SDEC,CDECZ,CDEC,
     :                 RADIF,SRADIF,CRADIF,DENOM

      DOUBLE PRECISION TINY
      PARAMETER (TINY=1D-6)


*  Trig functions
      SDECZ=SIN(DECZ)
      SDEC=SIN(DEC)
      CDECZ=COS(DECZ)
      CDEC=COS(DEC)
      RADIF=RA-RAZ
      SRADIF=SIN(RADIF)
      CRADIF=COS(RADIF)

*  Reciprocal of star vector length to tangent plane
      DENOM=SDEC*SDECZ+CDEC*CDECZ*CRADIF

*  Handle vectors too far from axis
      IF (DENOM.GT.TINY) THEN
         J=0
      ELSE IF (DENOM.GE.0D0) THEN
         J=1
         DENOM=TINY
      ELSE IF (DENOM.GT.-TINY) THEN
         J=2
         DENOM=-TINY
      ELSE
         J=3
      END IF

*  Compute tangent plane coordinates (even in dubious cases)
      XI=CDEC*SRADIF/DENOM
      ETA=(SDEC*CDECZ-CDEC*SDECZ*CRADIF)/DENOM

      END
