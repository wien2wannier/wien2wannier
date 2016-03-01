!!! wien2wannier/SRC_wplot/moduls.f

module param
  use const, only: BUFSZ, DPk, R8, C16, PI

  implicit none

  public

  character(*), parameter, private :: rev_str="$version: v1.0.0-97-g46c0317$"
  character(*), parameter, public  :: &
       wien2wannier_version = rev_str(11 : len (rev_str)-1)

!!!/=== The following is from the old ‘param.inc_{r,c}’ ===================
!     >> Standard WIEN97 Parameters <<
!
!     Constant parameter definition
!
  integer, parameter :: Nrad=981, Nsym=48, lmax7=8, LOmax=3, NRF=4, NLOat=3

!:17[
!     There are a couple of optional features in how to set-up the wave
!     function's augmentation coefficients. None of them are essential
!     for the wave functions themselves. However, to stay consistent with
!     LAPW1 and LAPW2 they have been added to LAPW7 as well. To allow
!     access to these optional features four logical control parameters
!     have been introduced which have to be set to .TRUE. to activate
!     the LAPW1/LAPW2 features. Any other choice of these parameters
!     is possible as well, but the code has not been optimized for any
!     of these alternative choices. For more details see the "trick 17"
!     sections C:17[ ... C:17] in the code (Uwe Birkenheuer)
!
  logical, parameter :: &
       kconjg=.true., UseRot=.true., AddLoc=.true., mvatom=.true.
!:17]
!!!!!\=== END param.inc ===================================================

  real(DPk), parameter :: clight=137.0359991_DPk

  integer, parameter :: unit_in=5, unit_out=6, unit_struct=8, unit_vector=10
  integer, parameter :: unit_grid=7, unit_vsp=18, unit_inwf=31, unit_chk=32
  integer, parameter :: unit_psink=21, unit_psiarg=22, unit_def=1
  integer, parameter :: unit_rot=33

!!! The following is set by ‘wplot.f’ to be passed to ‘main.F’
  character(BUFSZ) :: vecfn, psinkfn, psiargfn, outfn
  integer          :: idx_wann=0, iproc=0
end module param

module ams
  use const, only: DPk

  ! FIXME: get more precise atomic mass data
  real(DPk), parameter :: atom_mass(103) = (/ &
       &    1.0,   4.0,   6.9,   9.0,  10.8,  12.0,  14.0,  16.0,  19.0, &
       &   20.2,  23.0,  24.3,  27.0,  28.1,  31.0,  32.0,  35.4,  40.0, &
       &   39.1,  40.0,  45.0,  47.9,  50.9,  52.0,  54.9,  55.8,  58.9, &
       &   58.7,  63.5,  65.4,  69.7,  72.6,  74.9,  79.0,  79.9,  83.8, &
       &   85.5,  87.6,  88.9,  91.2,  92.9,  95.9,  98.0, 101.1, 102.9, &
       &  106.4, 107.9, 112.4, 114.8, 118.7, 121.8, 127.6, 126.9, 131.3, &
       &  132.9, 137.3, 138.9, 140.1, 140.9, 144.2, 145.0, 150.4, 152.0, &
       &  157.3, 158.9, 162.5, 164.9, 167.3, 168.9, 173.0, 175.0, 178.5, &
       &  180.9, 183.8, 186.2, 190.2, 192.2, 195.1, 197.0, 200.6, 204.4, &
       &  207.2, 209.0, 209.0, 210.0, 222.0, 223.0, 226.0, 227.0, 232.0, &
       &  231.0, 238.0, 237.0, 244.0, 243.0, 247.0, 247.0, 251.0, 252.0, &
       &  257.0, 258.0, 259.0, 262.0                                     &
       /)
end module ams


!!! “Minimodules” (converted COMMON blocks &c.)
module atspdt
  use const, only: DPk
  implicit none

  real(DPk),allocatable :: P(:,:,:),DP(:,:,:)
end module atspdt

module bessfu
  use const, only: DPk
  implicit none

  real(DPk),allocatable ::  FJ(:,:,:),DFJ(:,:,:),RAD(:)
  integer,  allocatable :: IRAD(:)
end module bessfu

module grid
  use const, only: DPk
  implicit none
  real(DPk),allocatable :: rgrid(:,:)
  integer,  allocatable :: ireg(:),ilat(:,:),iri(:)
  integer npg
end module grid

module latt
  use const, only: DPk

  implicit none

  private DPk
  public

! BR1(i,:) -- the real space lattice vectors a_i of the conventional u.c.
! BR2(i,:) -- the real space lattice vectors a_i of the primitive unit cell
! BR3(i,:) -- the reciprocal lattice vectors b_i of the conventional u.c
! BR4(i,:) -- the reciprocal lattice vectors b_i of the primitive unit cell


  real(DPk) :: BR1(3,3), BR2(3,3), BR3(3,3), BR4(3,3)
end module latt

module loabc
  use const, only: DPk
  implicit none

  real(DPk),allocatable :: ALO(:,:,:,:)
end module loabc

module lolog
  implicit none

  integer :: NLO
  logical, allocatable :: LAPW(:,:)
  integer, allocatable :: ILO(:,:)
end module lolog

module ps1
  use const, only: DPK
  implicit none

  real(DPk) :: dep(5),deq(5),dd,dvc,dsal,dk,dm1
end module ps1

module radfu
  use const, only: DPk
  implicit none

  real(DPk),allocatable ::  RRAD(:,:,:,:)
end module radfu

module radgrd
  use const, only: DPk
  implicit none

  real(DPk),allocatable :: RM(:,:), DX(:)
end module radgrd

module struct
  use const, only: DPk
  implicit none

  real(DPk), allocatable :: POS(:,:), RMT(:)
end module struct

module sym2
  use const, only: DPk
  use param, only: NSYM

  implicit none

  private :: NSYM
  public

  real(DPk) :: TAU (3,NSYM)
  integer   :: IMAT(3,3,NSYM), IORD
end module sym2

module work
  use const, only: DPK
  implicit none

  complex(DPK), allocatable :: aug(:,:,:)
end module work

module work1
  use const, only: DPK
  use param, only: Nrad
  implicit none

  real(DPk) :: A(Nrad), B(Nrad)
end module work1

!!/---
!! Local Variables:
!! mode: f90
!! End:
!!\---
!!
!! Time-stamp: <2016-03-01 11:35:12 assman@faepop36.tu-graz.ac.at>
