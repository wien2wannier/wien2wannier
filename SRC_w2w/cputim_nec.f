!!! wien2wannier/SRC_w2w/cputim_nec.f
!!!
!!! $Id: cputim_nec.f 167 2014-02-03 09:43:33Z assmann $

      SUBROUTINE CPUTIM(TIME)
      DOUBLE PRECISION   TIME
!
      CALL CLOCK(TIME)
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
