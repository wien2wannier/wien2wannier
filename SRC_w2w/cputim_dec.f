!!! wien2wannier/SRC_w2w/cputim_dec.f
!!!
!!! $Id: cputim_dec.f 385 2015-06-01 13:08:18Z assmann $

	subroutine cputim(cp)
	real*8 cp
	integer tms (4)
	integer utime,stime
	equivalence (tms(1),utime)
	equivalence (tms(2),stime)
	integer cutime,cstime
	equivalence (tms(3),cutime)
	equivalence (tms(4),cstime)

	integer HZ
	parameter (HZ = 100)

	call time (tms)
	cp = dble (secnds(0.0))
	return
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
