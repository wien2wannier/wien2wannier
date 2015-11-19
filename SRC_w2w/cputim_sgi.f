!!! wien2wannier/SRC_w2w/cputim_sgi.f
!!!
!!! $Id: cputim_sgi.f 385 2015-06-01 13:08:18Z assmann $

      SUBROUTINE CPUTIM(TIME)
      real*8 time                                                  
      time=second()
      end
      function second(dummy)
      real tarray(2)
      second=etime(tarray)
      end
      SUBROUTINE WALLTIM(DSEC)
      REAL*8 DSEC
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
