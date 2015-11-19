!!! wien2wannier/SRC_w2w/cputim_cray.f
!!!
!!! $Id: cputim_cray.f 385 2015-06-01 13:08:18Z assmann $

      SUBROUTINE CPUTIM(DSEC)
!
!        Scalar Arguments
!
      DOUBLE PRECISION   DSEC
!
!        Locals
!
      INTEGER            ISEC
      INTEGER            DUMMY(4)
!
!        External Functions
!
      INTEGER            TIMES
      EXTERNAL           TIMES
!
!        Intrinsic Functions
!
      INTRINSIC          DBLE
!
      ISEC = TIMES(DUMMY)
      DSEC = DBLE(DUMMY(1))/240000000.0D0
!
      RETURN
!
!        End of 'CPUTIM'
!
      END

      SUBROUTINE WALLTIM(DSEC)
      DOUBLE PRECISION DSEC
      DSEC=0.0D0
      RETURN
      END


!!/---
!! Local Variables:
!! mode: f90
!! End:
!!\---
!!
!! Time-stamp: <2015-05-23 19:58:48 elias>
