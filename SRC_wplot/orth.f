!!! wien2wannier/SRC_wplot/orth.f
!!!
!!! $Id: orth.f 385 2015-06-01 13:08:18Z assmann $

LOGICAL FUNCTION ORTH(A)
  use const
  IMPLICIT REAL(R8) (A-H,O-Z) 
  DIMENSION A(3,3)
  !
  ! Function to check whether A is orthogonal or not
  ! 
  DATA TOL / 1.0D-4 /

  ORTH = .FALSE.
  DO I=1,3
     DO J=1,3
        T = A(1,I)*A(1,J) + A(2,I)*A(2,J) + A(3,I)*A(3,J)
        IF(J.EQ.I) T = T - 1.0D0
        IF(ABS(T).GT.TOL) RETURN
     END DO
  END DO
  ORTH = .TRUE.
  RETURN
END FUNCTION ORTH


!!/---
!! Local Variables:
!! mode: f90
!! End:
!!\---
!!
!! Time-stamp: <2015-05-23 19:58:48 elias>
