!!! wien2wannier/SRC_w2w/mod_rc.F
!!!
!!!    Modules for wien2wannier.  This file contains modules that
!!!    depend on real/complex compilation.
!!!
!!! Copyright 2010-2012 Jan Kuneš, Philipp Wissgott
!!!           2013-2016 Elias Assmann

!!/=== Defined here: =============================
!!
!! xa3: init_xa3(), BKX(:),BKY(:),BKZ(:), GX(:,:),GY(:,:),GZ(:,:),
!!      XK(:),YK(:),ZK(:), vecsz(:)
!!
!!\===============================================

module xa3
  use const, only: R8
#ifdef _COMPLEX_
  use const, only: C16
#endif

  implicit none
  private
  public :: init_xa3, XK,YK,ZK, BKX,BKY,BKZ, GX,GY,GZ, vecsz, a

#ifndef _COMPLEX_
  real(R8), &
#else
  complex(C16), &
#endif
  allocatable :: a(:,:,:)

  integer,  allocatable :: GX(:,:),GY(:,:),GZ(:,:), vecsz(:)
  real(R8), allocatable :: BKX(:), BKY(:), BKZ(:)
  real(R8), allocatable :: XK(:),  YK(:),  ZK(:)

contains
  subroutine init_xa3(nb, nmat, num_kpts)
    integer, intent(in) :: NB, Nmat, num_kpts

    allocate( a  (Nmat, Nb, num_kpts), vecsz(num_kpts) )
    allocate( BKX(Nmat),         BKY(Nmat),         BKZ(Nmat),          &
         &    GX (Nmat,num_kpts), GY(Nmat,num_kpts), GZ(Nmat,num_kpts), &
         &    XK (     num_kpts), YK(     num_kpts), ZK(     num_kpts))

    GX = 0; GY = 0; GZ = 0
  end subroutine init_xa3
end module xa3

!!/---
!! Local Variables:
!! mode: f90
!! End:
!!\---
!!
!! Time-stamp: <2016-07-18 17:07:11 assman@faepop71.tu-graz.ac.at>
