!!! wien2wannier/SRC_wplot/wavint.f

SUBROUTINE WAVINT(R,NPW,PSI,bk,coef,nmat)
  use const, only: DPk

  implicit none

  integer,      intent(in)  :: NPW, NMat
  real(DPk),    intent(in)  :: R(3), BK(3,Nmat)
  complex(DPk), intent(out) :: Psi
  complex(DPk), intent(in)  :: coef(NMat)

! evaluation of the wave function in the interstitial region:
! psi(r) = Sum(K) c_K/sqrt(V) e^i(K+k)r
! --------------------------------------------------------------
! Input:
! R    -- grid point in (gobal) Cartesian coordinates
! NPW  -- number of PW basis functions
!
! MODULE EIGVEC
! BK   -- the PW wave vectors K+k in (gobal) Cartesian coordinates
! COEF -- the PW coefficients c_K (including 1/sqrt(V))
!
! Output:
! PSI  -- the wavr function psi(r)
! --------------------------------------------------------------

  real(DPk) :: arg
  integer   :: iPW

  PSI = 0
  DO IPW=1,NPW
     ARG = R(1)*BK(1,IPW) + R(2)*BK(2,IPW) + R(3)*BK(3,IPW)
     PSI = PSI + COEF(IPW) * cmplx(cos(arg), sin(arg), DPk)
  END DO
END SUBROUTINE WAVINT


!!/---
!! Local Variables:
!! mode: f90
!! End:
!!\---
!!
!! Time-stamp: <2016-07-15 11:44:47 assman@faepop71.tu-graz.ac.at>
