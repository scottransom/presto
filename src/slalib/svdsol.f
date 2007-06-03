      SUBROUTINE sla_SVDSOL (M, N, MP, NP, B, U, W, V, WORK, X)
*+
*     - - - - - - -
*      S V D S O L
*     - - - - - - -
*
*  From a given vector and the SVD of a matrix (as obtained from
*  the SVD routine), obtain the solution vector (double precision)
*
*  This routine solves the equation:
*
*     A . x = b
*
*  where:
*
*     A   is a given M (rows) x N (columns) matrix, where M.GE.N
*     x   is the N-vector we wish to find
*     b   is a given M-vector
*
*  by means of the Singular Value Decomposition method (SVD).  In
*  this method, the matrix A is first factorised (for example by
*  the routine sla_SVD) into the following components:
*
*     A = U x W x VT
*
*  where:
*
*     A   is the M (rows) x N (columns) matrix
*     U   is an M x N column-orthogonal matrix
*     W   is an N x N diagonal matrix with W(I,I).GE.0
*     VT  is the transpose of an NxN orthogonal matrix
*
*     Note that M and N, above, are the LOGICAL dimensions of the
*     matrices and vectors concerned, which can be located in
*     arrays of larger PHYSICAL dimensions MP and NP.
*
*  The solution is found from the expression:
*
*     x = V . [diag(1/Wj)] . (transpose(U) . b)
*
*  Notes:
*
*  1)  If matrix A is square, and if the diagonal matrix W is not
*      adjusted, the method is equivalent to conventional solution
*      of simultaneous equations.
*
*  2)  If M>N, the result is a least-squares fit.
*
*  3)  If the solution is poorly determined, this shows up in the
*      SVD factorisation as very small or zero Wj values.  Where
*      a Wj value is small but non-zero it can be set to zero to
*      avoid ill effects.  The present routine detects such zero
*      Wj values and produces a sensible solution, with highly
*      correlated terms kept under control rather than being allowed
*      to elope to infinity, and with meaningful values for the
*      other terms.
*
*  Given:
*     M,N    i         numbers of rows and columns in matrix A
*     MP,NP  i         physical dimensions of array containing matrix A
*     B      d(M)      known vector b
*     U      d(MP,NP)  array containing MxN matrix U
*     W      d(N)      NxN diagonal matrix W (diagonal elements only)
*     V      d(NP,NP)  array containing NxN orthogonal matrix V
*
*  Returned:
*     WORK   d(N)      workspace
*     X      d(N)      unknown vector x
*
*  Reference:
*     Numerical Recipes, section 2.9.
*
*  P.T.Wallace   Starlink   29 October 1993
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

      INTEGER M,N,MP,NP
      DOUBLE PRECISION B(M),U(MP,NP),W(N),V(NP,NP),WORK(N),X(N)

      INTEGER J,I,JJ
      DOUBLE PRECISION S



*  Calculate [diag(1/Wj)] . transpose(U) . b (or zero for zero Wj)
      DO J=1,N
         S=0D0
         IF (W(J).NE.0D0) THEN
            DO I=1,M
               S=S+U(I,J)*B(I)
            END DO
            S=S/W(J)
         END IF
         WORK(J)=S
      END DO

*  Multiply by matrix V to get result
      DO J=1,N
         S=0D0
         DO JJ=1,N
            S=S+V(J,JJ)*WORK(JJ)
         END DO
         X(J)=S
      END DO

      END
