!!! wien2wannier/SRC_w2w/ph.f

      REAL(R8) FUNCTION PH(N)
      use const
      IMPLICIT REAL(R8) (A-H,O-Z)
      PH=1.0D0
      K=MOD(IABS(N),4)
      IF(K.EQ.0) GOTO 5
      PH=-1.0D0
    5 RETURN
      END


!!/---
!! Local Variables:
!! mode: f90
!! End:
!!\---
!!
!! Time-stamp: <2015-05-23 19:58:48 elias>
