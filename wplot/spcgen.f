      SUBROUTINE SPCGEN(NAT,RMT,ATMS)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION RMT(NAT),ATMS(3,NAT)
!
! -----------------------------------------------------------------------
! For a given lattice basis {a_1,a_2,a_3} each muffin tin sphere 
!   MTS := { r | |r| < Rmt }
! is surrounded by a (smallest) parallel epiped
!   SPC := Sum(i=1,3) [-s_i,+s_i] * a_i
! with s_i > 0 for all three (primitive) lattice vectors a_i  .
!
! The mininal realization of such a parallel epiped is given by 
!   s_i = Rmt * |b_i|
! where {b_1,b_2,b_3} is the reciprocal lattice basis of {a_1,a_2,a_3}.
! -----------------------------------------------------------------------
!
! Input:
! NAT       -- number of inequivalent atoms
! RMT(j)    -- muffin tin radius of the j-th inequivalent atom
!
! Output:
! ATMS(:,j) -- the size parameter s_i of the j-th inequivalent atom
!
      COMMON /LATT/ VUC,BR1(3,3),BR2(3,3),BR3(3,3),BR4(3,3)
!
      DIMENSION B(3)
!
      DO 10 I=1,3
        B(I) = SQRT( BR4(I,1)**2 + BR4(I,2)**2 + BR4(I,3)**2 )
   10 CONTINUE
      DO 20 JATOM=1,NAT
        DO 20 I=1,3
          ATMS(I,JATOM) = RMT(JATOM) * B(I)
   30   CONTINUE
   20 CONTINUE
!
      END
