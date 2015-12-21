!!! wien2wannier/SRC_wplot/spcgen.f

SUBROUTINE SPCGEN(NAT,RMT,ATMS)
  use const
  use latt
  IMPLICIT REAL(R8) (A-H,O-Z)
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

  DIMENSION B(3)

  DO I=1,3
     B(I) = SQRT( BR4(I,1)**2 + BR4(I,2)**2 + BR4(I,3)**2 )
  END DO
  DO JATOM=1,NAT
     DO I=1,3
        ATMS(I,JATOM) = RMT(JATOM) * B(I)
     END DO
  END DO

END SUBROUTINE SPCGEN


!!/---
!! Local Variables:
!! mode: f90
!! End:
!!\---
!!
!! Time-stamp: <2015-05-23 19:58:48 elias>
