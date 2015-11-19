      LOGICAL FUNCTION ORTH(A)
      IMPLICIT REAL*8 (A-H,O-Z) 
      DIMENSION A(3,3)
!
! Function to check whether A is orthogonal or not
! 
      DATA TOL / 1.0D-4 /

      ORTH = .FALSE.
      DO 10 I=1,3
        DO 20 J=1,3
          T = A(1,I)*A(1,J) + A(2,I)*A(2,J) + A(3,I)*A(3,J)
          IF(J.EQ.I) T = T - 1.0D0
          IF(ABS(T).GT.TOL) RETURN
   20   CONTINUE
   10 CONTINUE
      ORTH = .TRUE.
      RETURN
      END
 
       
