!!! wien2wannier/SRC_w2w/rout.f
!!!
!!! $Id: rout.f 385 2015-06-01 13:08:18Z assmann $

      SUBROUTINE ROUT(LL)
      USE param

      use const, only: R8

      IMPLICIT REAL(R8) (A-H,O-Z)

      COMMON /RINTEG/  RI_MAT(0:lmax2,nrf,nrf,2,2)

      write(unit_out,*)'LL=',ll
      WRITE(unit_out,*)'RADIAL INTEGRALS:'
      DO 10 IRF1=1,NRF
      DO 11 IS1=1,2
      DO 12 IS2=1,2
      DO 13 IRF2=1,NRF
      if (abs(RI_MAT(LL,IRF1,IRF2,IS1,IS2)).gt.1d-4) &
      WRITE(unit_out,*)irf1,irf2,'     ',is1,is2,RI_MAT(LL,IRF1,IRF2,IS1,IS2)
 13   continue 
 12   continue 
 11   continue 
 10   continue 
      end


!!/---
!! Local Variables:
!! mode: f90
!! End:
!!\---
!!
!! Time-stamp: <2015-05-23 19:58:48 elias>
