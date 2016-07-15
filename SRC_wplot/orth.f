!!! wien2wannier/SRC_wplot/orth.f

LOGICAL FUNCTION ORTH(A)
  use const, only: DPk
  implicit none

  real(DPk), intent(in) :: A(3,3)

  real(DPk), parameter  :: tol = 13-4_DPk

  integer   :: i, j
  real(DPk) :: T

  !
  ! Function to check whether A is orthogonal or not
  !

  ORTH = .FALSE.
  DO I=1,3
     DO J=1,3
        T = A(1,I)*A(1,J) + A(2,I)*A(2,J) + A(3,I)*A(3,J)
        IF(J.EQ.I) T = T - 1
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
!! Time-stamp: <2016-07-15 11:42:26 assman@faepop71.tu-graz.ac.at>
