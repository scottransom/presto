      SUBROUTINE sla_SVDCOV (N, NP, NC, W, V, WORK, CVM)
*+
*     - - - - - - -
*      S V D C O V
*     - - - - - - -
*
*  From the W and V matrices from the SVD factorisation of a matrix
*  (as obtained from the sla_SVD routine), obtain the covariance matrix.
*
*  (double precision)
*
*  Given:
*     N      i         number of rows and columns in matrices W and V
*     NP     i         first dimension of array containing matrix V
*     NC     i         first dimension of array to receive CVM
*     W      d(N)      NxN diagonal matrix W (diagonal elements only)
*     V      d(NP,NP)  array containing NxN orthogonal matrix V
*
*  Returned:
*     WORK   d(N)      workspace
*     CVM    d(NC,NC)  array to receive covariance matrix
*
*  Reference:
*     Numerical Recipes, section 14.3.
*
*  P.T.Wallace   Starlink   December 1988
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

      INTEGER N,NP,NC
      DOUBLE PRECISION W(N),V(NP,NP),WORK(N),CVM(NC,NC)

      INTEGER I,J,K
      DOUBLE PRECISION S



      DO I=1,N
         S=W(I)
         IF (S.NE.0D0) THEN
            WORK(I)=1D0/(S*S)
         ELSE
            WORK(I)=0D0
         END IF
      END DO
      DO I=1,N
         DO J=1,I
            S=0D0
            DO K=1,N
               S=S+V(I,K)*V(J,K)*WORK(K)
            END DO
            CVM(I,J)=S
            CVM(J,I)=S
         END DO
      END DO

      END
