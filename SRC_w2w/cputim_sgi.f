!!! wien2wannier/SRC_w2w/cputim_sgi.f

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
