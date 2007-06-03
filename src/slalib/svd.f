      SUBROUTINE sla_SVD (M, N, MP, NP, A, W, V, WORK, JSTAT)
*+
*     - - - -
*      S V D
*     - - - -
*
*  Singular value decomposition  (double precision)
*
*  This routine expresses a given matrix A as the product of
*  three matrices U, W, V:
*
*     A = U x W x VT
*
*  Where:
*
*     A   is any M (rows) x N (columns) matrix, where M.GE.N
*     U   is an M x N column-orthogonal matrix
*     W   is an N x N diagonal matrix with W(I,I).GE.0
*     VT  is the transpose of an N x N orthogonal matrix
*
*     Note that M and N, above, are the LOGICAL dimensions of the
*     matrices and vectors concerned, which can be located in
*     arrays of larger PHYSICAL dimensions, given by MP and NP.
*
*  Given:
*     M,N    i         numbers of rows and columns in matrix A
*     MP,NP  i         physical dimensions of array containing matrix A
*     A      d(MP,NP)  array containing MxN matrix A
*
*  Returned:
*     A      d(MP,NP)  array containing MxN column-orthogonal matrix U
*     W      d(N)      NxN diagonal matrix W (diagonal elements only)
*     V      d(NP,NP)  array containing NxN orthogonal matrix V
*     WORK   d(N)      workspace
*     JSTAT  i         0 = OK, -1 = A wrong shape, >0 = index of W
*                      for which convergence failed.  See note 2, below.
*
*   Notes:
*
*   1)  V contains matrix V, not the transpose of matrix V.
*
*   2)  If the status JSTAT is greater than zero, this need not
*       necessarily be treated as a failure.  It means that, due to
*       chance properties of the matrix A, the QR transformation
*       phase of the routine did not fully converge in a predefined
*       number of iterations, something that very seldom occurs.
*       When this condition does arise, it is possible that the
*       elements of the diagonal matrix W have not been correctly
*       found.  However, in practice the results are likely to
*       be trustworthy.  Applications should report the condition
*       as a warning, but then proceed normally.
*
*  References:
*     The algorithm is an adaptation of the routine SVD in the EISPACK
*     library (Garbow et al 1977, EISPACK Guide Extension, Springer
*     Verlag), which is a FORTRAN 66 implementation of the Algol
*     routine SVD of Wilkinson & Reinsch 1971 (Handbook for Automatic
*     Computation, vol 2, ed Bauer et al, Springer Verlag).  These
*     references give full details of the algorithm used here.  A good
*     account of the use of SVD in least squares problems is given in
*     Numerical Recipes (Press et al 1986, Cambridge University Press),
*     which includes another variant of the EISPACK code.
*
*  Last revision:   8 September 2005
*
*  Copyright P.T.Wallace.  All rights reserved.
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
      DOUBLE PRECISION A(MP,NP),W(N),V(NP,NP),WORK(N)
      INTEGER JSTAT

*  Maximum number of iterations in QR phase
      INTEGER ITMAX
      PARAMETER (ITMAX=30)

      INTEGER L,L1,I,K,J,K1,ITS,I1
      LOGICAL CANCEL
      DOUBLE PRECISION G,SCALE,AN,S,X,F,H,C,Y,Z



*  Variable initializations to avoid compiler warnings.
      L = 0
      L1 = 0

*  Check that the matrix is the right shape
      IF (M.LT.N) THEN

*     No:  error status
         JSTAT = -1

      ELSE

*     Yes:  preset the status to OK
         JSTAT = 0

*
*     Householder reduction to bidiagonal form
*     ----------------------------------------

         G = 0D0
         SCALE = 0D0
         AN = 0D0
         DO I=1,N
            L = I+1
            WORK(I) = SCALE*G
            G = 0D0
            S = 0D0
            SCALE = 0D0
            IF (I.LE.M) THEN
               DO K=I,M
                  SCALE = SCALE+ABS(A(K,I))
               END DO
               IF (SCALE.NE.0D0) THEN
                  DO K=I,M
                     X = A(K,I)/SCALE
                     A(K,I) = X
                     S = S+X*X
                  END DO
                  F = A(I,I)
                  G = -SIGN(SQRT(S),F)
                  H = F*G-S
                  A(I,I) = F-G
                  IF (I.NE.N) THEN
                     DO J=L,N
                        S = 0D0
                        DO K=I,M
                           S = S+A(K,I)*A(K,J)
                        END DO
                        F = S/H
                        DO K=I,M
                           A(K,J) = A(K,J)+F*A(K,I)
                        END DO
                     END DO
                  END IF
                  DO K=I,M
                     A(K,I) = SCALE*A(K,I)
                  END DO
               END IF
            END IF
            W(I) = SCALE*G
            G = 0D0
            S = 0D0
            SCALE = 0D0
            IF (I.LE.M .AND. I.NE.N) THEN
               DO K=L,N
                  SCALE = SCALE+ABS(A(I,K))
               END DO
               IF (SCALE.NE.0D0) THEN
                  DO K=L,N
                     X = A(I,K)/SCALE
                     A(I,K) = X
                     S = S+X*X
                  END DO
                  F = A(I,L)
                  G = -SIGN(SQRT(S),F)
                  H = F*G-S
                  A(I,L) = F-G
                  DO K=L,N
                     WORK(K) = A(I,K)/H
                  END DO
                  IF (I.NE.M) THEN
                     DO J=L,M
                        S = 0D0
                        DO K=L,N
                           S = S+A(J,K)*A(I,K)
                        END DO
                        DO K=L,N
                           A(J,K) = A(J,K)+S*WORK(K)
                        END DO
                     END DO
                  END IF
                  DO K=L,N
                     A(I,K) = SCALE*A(I,K)
                  END DO
               END IF
            END IF

*        Overestimate of largest column norm for convergence test
            AN = MAX(AN,ABS(W(I))+ABS(WORK(I)))

         END DO

*
*     Accumulation of right-hand transformations
*     ------------------------------------------

         DO I=N,1,-1
            IF (I.NE.N) THEN
               IF (G.NE.0D0) THEN
                  DO J=L,N
                     V(J,I) = (A(I,J)/A(I,L))/G
                  END DO
                  DO J=L,N
                     S = 0D0
                     DO K=L,N
                        S = S+A(I,K)*V(K,J)
                     END DO
                     DO K=L,N
                        V(K,J) = V(K,J)+S*V(K,I)
                     END DO
                  END DO
               END IF
               DO J=L,N
                  V(I,J) = 0D0
                  V(J,I) = 0D0
               END DO
            END IF
            V(I,I) = 1D0
            G = WORK(I)
            L = I
         END DO

*
*     Accumulation of left-hand transformations
*     -----------------------------------------

         DO I=N,1,-1
            L = I+1
            G = W(I)
            IF (I.NE.N) THEN
               DO J=L,N
                  A(I,J) = 0D0
               END DO
            END IF
            IF (G.NE.0D0) THEN
               IF (I.NE.N) THEN
                  DO J=L,N
                     S = 0D0
                     DO K=L,M
                        S = S+A(K,I)*A(K,J)
                     END DO
                     F = (S/A(I,I))/G
                     DO K=I,M
                        A(K,J) = A(K,J)+F*A(K,I)
                     END DO
                  END DO
               END IF
               DO J=I,M
                  A(J,I) = A(J,I)/G
               END DO
            ELSE
               DO J=I,M
                  A(J,I) = 0D0
               END DO
            END IF
            A(I,I) = A(I,I)+1D0
         END DO

*
*     Diagonalisation of the bidiagonal form
*     --------------------------------------

         DO K=N,1,-1
            K1 = K-1

*        Iterate until converged
            ITS = 0
            DO WHILE (ITS.LT.ITMAX)
               ITS = ITS+1

*           Test for splitting into submatrices
               CANCEL = .TRUE.
               DO L=K,1,-1
                  L1 = L-1
                  IF (AN+ABS(WORK(L)).EQ.AN) THEN
                     CANCEL = .FALSE.
                     GO TO 10
                  END IF
*              (Following never attempted for L=1 because WORK(1) is zero)
                  IF (AN+ABS(W(L1)).EQ.AN) GO TO 10
               END DO
 10            CONTINUE

*           Cancellation of WORK(L) if L>1
               IF (CANCEL) THEN
                  C = 0D0
                  S = 1D0
                  DO I=L,K
                     F = S*WORK(I)
                     IF (AN+ABS(F).EQ.AN) GO TO 20
                     G = W(I)
                     H = SQRT(F*F+G*G)
                     W(I) = H
                     C = G/H
                     S = -F/H
                     DO J=1,M
                        Y = A(J,L1)
                        Z = A(J,I)
                        A(J,L1) = Y*C+Z*S
                        A(J,I) = -Y*S+Z*C
                     END DO
                  END DO
 20               CONTINUE
               END IF

*           Converged?
               Z = W(K)
               IF (L.EQ.K) THEN

*              Yes:  stop iterating
                  ITS = ITMAX

*              Ensure singular values non-negative
                  IF (Z.LT.0D0) THEN
                     W(K) = -Z
                     DO J=1,N
                        V(J,K) = -V(J,K)
                     END DO
                  END IF
               ELSE

*              Not converged yet:  set status if iteration limit reached
                  IF (ITS.EQ.ITMAX) JSTAT = K

*              Shift from bottom 2x2 minor
                  X = W(L)
                  Y = W(K1)
                  G = WORK(K1)
                  H = WORK(K)
                  F = ((Y-Z)*(Y+Z)+(G-H)*(G+H))/(2D0*H*Y)
                  IF (ABS(F).LE.1D15) THEN
                     G = SQRT(F*F+1D0)
                  ELSE
                     G = ABS(F)
                  END IF
                  F = ((X-Z)*(X+Z)+H*(Y/(F+SIGN(G,F))-H))/X

*              Next QR transformation
                  C = 1D0
                  S = 1D0
                  DO I1=L,K1
                     I = I1+1
                     G = WORK(I)
                     Y = W(I)
                     H = S*G
                     G = C*G
                     Z = SQRT(F*F+H*H)
                     WORK(I1) = Z
                     IF (Z.NE.0D0) THEN
                        C = F/Z
                        S = H/Z
                     ELSE
                        C = 1D0
                        S = 0D0
                     END IF
                     F = X*C+G*S
                     G = -X*S+G*C
                     H = Y*S
                     Y = Y*C
                     DO J=1,N
                        X = V(J,I1)
                        Z = V(J,I)
                        V(J,I1) = X*C+Z*S
                        V(J,I) = -X*S+Z*C
                     END DO
                     Z = SQRT(F*F+H*H)
                     W(I1) = Z
                     IF (Z.NE.0D0) THEN
                        C = F/Z
                        S = H/Z
                     END IF
                     F = C*G+S*Y
                     X = -S*G+C*Y
                     DO J=1,M
                        Y = A(J,I1)
                        Z = A(J,I)
                        A(J,I1) = Y*C+Z*S
                        A(J,I) = -Y*S+Z*C
                     END DO
                  END DO
                  WORK(L) = 0D0
                  WORK(K) = F
                  W(K) = X
               END IF
            END DO
         END DO
      END IF

      END
