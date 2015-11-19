!!! wien2wannier/SRC_w2w/cputim_sgi.f
!!!
!!! $Id: cputim_sgi.f 167 2014-02-03 09:43:33Z assmann $

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
