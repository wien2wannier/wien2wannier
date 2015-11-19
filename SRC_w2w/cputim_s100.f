!!! wien2wannier/SRC_w2w/cputim_s100.f
!!!
!!! $Id: cputim_s100.f 385 2015-06-01 13:08:18Z assmann $

      SUBROUTINE CPUTIM(TIME)
      DOUBLE PRECISION   TIME
!
      CALL CLOCK(TIME,0,2)
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
