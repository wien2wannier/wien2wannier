!!! wien2wannier/SRC_wplot/wavint.f

SUBROUTINE WAVINT(R,NPW,PSI,bk,coef,nmat)
  use const
  use param

  IMPLICIT REAL(R8) (A-H,O-Z)
  DIMENSION  R(3)
  COMPLEX(C16) PSI
!
! evaluation of the wave function in the interstitial region:
! psi(r) = Sum(K) c_K/sqrt(V) e^i(K+k)r 
! --------------------------------------------------------------
! Input:
! R    -- grid point in (gobal) Cartesian coordinates
! NPW  -- number of PW basis functions
!
! COMMON /EIGVEC/
! BK   -- the PW wave vectors K+k in (gobal) Cartesian coordinates
! COEF -- the PW coefficients c_K (including 1/sqrt(V))
!
! Output:
! PSI  -- the wavr function psi(r)
! --------------------------------------------------------------
  COMPLEX(C16) COEF(nmat) !changed by pwissgott

  real(r8)  BK(3,NMAT)

  PSI = (0.0D0,0.0D0)
  DO IPW=1,NPW
     ARG = R(1)*BK(1,IPW) + R(2)*BK(2,IPW) + R(3)*BK(3,IPW)
     PSI = PSI + COEF(IPW) * DCMPLX(COS(ARG),SIN(ARG))   !changed by P.wissgott
  END DO
END SUBROUTINE WAVINT


!!/---
!! Local Variables:
!! mode: f90
!! End:
!!\---
!!
!! Time-stamp: <2015-05-23 19:58:48 elias>
