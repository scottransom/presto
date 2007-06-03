      SUBROUTINE sla_SMAT (N, A, Y, D, JF, IW)
*+
*     - - - - -
*      S M A T
*     - - - - -
*
*  Matrix inversion & solution of simultaneous equations
*  (single precision)
*
*  For the set of n simultaneous equations in n unknowns:
*     A.Y = X
*
*  where:
*     A is a non-singular N x N matrix
*     Y is the vector of N unknowns
*     X is the known vector
*
*  SMATRX computes:
*     the inverse of matrix A
*     the determinant of matrix A
*     the vector of N unknowns
*
*  Arguments:
*
*     symbol  type dimension           before              after
*
*       N      int                 no. of unknowns       unchanged
*       A      real  (N,N)             matrix             inverse
*       Y      real   (N)              vector            solution
*       D      real                       -             determinant
*     * JF     int                        -           singularity flag
*       IW     int    (N)                 -              workspace
*
*  *  JF is the singularity flag.  If the matrix is non-singular,
*    JF=0 is returned.  If the matrix is singular, JF=-1 & D=0.0 are
*    returned.  In the latter case, the contents of array A on return
*    are undefined.
*
*  Algorithm:
*     Gaussian elimination with partial pivoting.
*
*  Speed:
*     Very fast.
*
*  Accuracy:
*     Fairly accurate - errors 1 to 4 times those of routines optimised
*     for accuracy.
*
*  Note:  replaces the obsolete sla_SMATRX routine.
*
*  P.T.Wallace   Starlink   10 September 1990
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

      INTEGER N
      REAL A(N,N),Y(N),D
      INTEGER JF
      INTEGER IW(N)

      REAL SFA
      PARAMETER (SFA=1E-20)

      INTEGER K,IMX,I,J,NP1MK,KI
      REAL AMX,T,AKK,YK,AIK


      JF=0
      D=1.0
      DO K=1,N
         AMX=ABS(A(K,K))
         IMX=K
         IF (K.NE.N) THEN
            DO I=K+1,N
               T=ABS(A(I,K))
               IF (T.GT.AMX) THEN
                  AMX=T
                  IMX=I
               END IF
            END DO
         END IF
         IF (AMX.LT.SFA) THEN
            JF=-1
         ELSE
            IF (IMX.NE.K) THEN
               DO J=1,N
                  T=A(K,J)
                  A(K,J)=A(IMX,J)
                  A(IMX,J)=T
               END DO
               T=Y(K)
               Y(K)=Y(IMX)
               Y(IMX)=T
               D=-D
            END IF
            IW(K)=IMX
            AKK=A(K,K)
            D=D*AKK
            IF (ABS(D).LT.SFA) THEN
               JF=-1
            ELSE
               AKK=1.0/AKK
               A(K,K)=AKK
               DO J=1,N
                  IF (J.NE.K) A(K,J)=A(K,J)*AKK
               END DO
               YK=Y(K)*AKK
               Y(K)=YK
               DO I=1,N
                  AIK=A(I,K)
                  IF (I.NE.K) THEN
                     DO J=1,N
                        IF (J.NE.K) A(I,J)=A(I,J)-AIK*A(K,J)
                     END DO
                     Y(I)=Y(I)-AIK*YK
                  END IF
               END DO
               DO I=1,N
                  IF (I.NE.K) A(I,K)=-A(I,K)*AKK
               END DO
            END IF
         END IF
      END DO
      IF (JF.NE.0) THEN
         D=0.0
      ELSE
         DO K=1,N
            NP1MK=N+1-K
            KI=IW(NP1MK)
            IF (NP1MK.NE.KI) THEN
               DO I=1,N
                  T=A(I,NP1MK)
                  A(I,NP1MK)=A(I,KI)
                  A(I,KI)=T
               END DO
            END IF
         END DO
      END IF
      END
