!!! wien2wannier/SRC_w2w/cputim_macg4.f
!!!
!!! $Id: cputim_macg4.f 167 2014-02-03 09:43:33Z assmann $

   subroutine cputim_(cp)
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
       IMPLICIT NONE
       DOUBLE PRECISION DSEC
       DSEC=0.0D0
       RETURN
       END

       SUBROUTINE WALLTIM_(DSEC)
       IMPLICIT NONE
       DOUBLE PRECISION DSEC
       DSEC=0.0D0
       RETURN
       END

!!/---
!! Local Variables:
!! mode: f90
!! End:
!!\---
