!!! wien2wannier/SRC_wplot/gbass.f
!!!
!!! $Id: gbass.f 166 2014-02-03 09:39:24Z assmann $

SUBROUTINE GBASS(RBAS,GBAS,TWOPI)
  use const
  IMPLICIT  REAL(R8) (A-H,O-Z)
  DIMENSION RBAS(3,3), GBAS(3,3)
  LOGICAL   TWOPI
  !
  ! << Input >>
  ! RBAS(i,:) -- the real space lattice vectros a_i
  ! TWOPI     -- normalization of the reciprocal lattice vectors:
  !                if .TRUE.  let < a_i | b_j > = 2pi * delta(i,j)
  !                if .FALSE. let < a_i | b_j > = delta(i,j)
  !
  ! << Output >>
  ! GBAS(j,:) -- the corresponding reciprocal lattice vectors b_j
  !
  !     << b_1 = a_2 x a_3 >>
  GBAS(1,1) = RBAS(2,2)*RBAS(3,3) - RBAS(2,3)*RBAS(3,2)
  GBAS(1,2) = RBAS(2,3)*RBAS(3,1) - RBAS(2,1)*RBAS(3,3)
  GBAS(1,3) = RBAS(2,1)*RBAS(3,2) - RBAS(2,2)*RBAS(3,1)
  !     << b_2 = a_3 x a_1 >>
  GBAS(2,1) = RBAS(3,2)*RBAS(1,3) - RBAS(3,3)*RBAS(1,2)
  GBAS(2,2) = RBAS(3,3)*RBAS(1,1) - RBAS(3,1)*RBAS(1,3)
  GBAS(2,3) = RBAS(3,1)*RBAS(1,2) - RBAS(3,2)*RBAS(1,1)
  !     << b_3 = a_1 x a_2 >>
  GBAS(3,1) = RBAS(1,2)*RBAS(2,3) - RBAS(1,3)*RBAS(2,2)
  GBAS(3,2) = RBAS(1,3)*RBAS(2,1) - RBAS(1,1)*RBAS(2,3)
  GBAS(3,3) = RBAS(1,1)*RBAS(2,2) - RBAS(1,2)*RBAS(2,1)
  !
  !     << normalization : VUC = < a_1 | a_2 x a_3 > = < a_1 | b_1 > >>
  VUC=0.0D0
  DO I=1,3
     VUC = VUC + RBAS(1,I)*GBAS(1,I)
  END DO

  IF(TWOPI)THEN
     FAC = PI / VUC
  ELSE
     FAC = 1  / VUC
  ENDIF

  GBAS = GBAS * FAC

END SUBROUTINE GBASS

!!/---
!! Local Variables:
!! mode: f90
!! End:
!!\---
