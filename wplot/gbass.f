      SUBROUTINE GBASS(RBAS,GBAS,TWOPI)
      IMPLICIT  REAL*8 (A-H,O-Z)
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
      DO 10 I=1,3
   10 VUC = VUC + RBAS(1,I)*GBAS(1,I)
      IF(TWOPI)THEN
        FAC = 3.14159265358979324D0 / VUC
      ELSE
        FAC = 1.0D0 / VUC
      ENDIF
      DO 20 I=1,3
      DO 20 J=1,3
   20 GBAS(I,J) = GBAS(I,J) * FAC
!
      RETURN
      END
