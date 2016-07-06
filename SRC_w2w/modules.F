!!! wien2wannier/SRC_w2w/modules.F

module param
  use const, only: R8, C16
  implicit none
  private :: R8, C16

  public

  character(*), parameter, private :: rev_str="$version: v1.0.0-122-gb2e6a61$"
  character(*), parameter, public  :: &
       wien2wannier_version = rev_str(11 : len (rev_str)-1)

  integer, parameter :: unit_in=5, unit_out=6, unit_amn=7, unit_fermi=51
  integer, parameter :: unit_vector=9, unit_nnkp=10, unit_eig=12, unit_vsp=18
  integer, parameter :: unit_struct=20, unit_ene=50, unit_def=1, unit_mmn=8

  !.....Optimize IBLOCK for your hardware (32-255)
  ! for ncom parameter check format 1003 in l2main.frc
  ! for x-dos set lxdos to 3
  integer, parameter :: IBLOCK=128, LMAX2=5, LOMAX=3, NCOM=49, NGAU=1500
  integer, parameter :: lxdos=3, Nloat=3, Nrad=881, Nrf=4

  integer :: NMAT=0, NDIF=0

  real(R8), parameter :: clight = 137.0359895_R8
end module param

module assleg
  use const, only: R8

  implicit none
  private

  public :: YR, N, init_assleg

  real(R8), allocatable :: YR(:,:)
  integer,  parameter   :: N=6, MAXDIM=81
contains
  subroutine init_assleg
    implicit none

    allocate( YR(N, MAXDIM) )
    YR=0
  end subroutine init_assleg
end module assleg

MODULE struct
  use const, only: R8
  implicit none
  private :: R8

  public

  logical                    :: rel
  real(R8)                   :: AA,BB,CC,VOL,pia(3),alpha(3)
  real(R8),allocatable       :: R0(:),DX(:),RMT(:),zz(:),rotloc(:,:,:),v(:)
  real(R8),allocatable       :: tau(:,:)
  real(R8),allocatable       :: pos(:,:)
  real(R8),allocatable       :: rotij(:,:,:),tauij(:,:)
  character(4)               :: lattic,irel,cform
  character(80)              :: title
  character(10), allocatable :: aname(:)
  integer                    :: nat,iord
  integer,allocatable        :: mult(:),jrj(:)
  integer,allocatable        :: iz(:,:,:)

 CONTAINS
  SUBROUTINE init_struct
    use param,      only: unit_struct, unit_out, Ndif
    use reallocate, only: realloc
    use const,      only: R8

    implicit none

    integer  :: ios, inum, isplit, iatnr
    real(R8) :: test
    integer  :: index, i, j, j1, j2, m, jatom

    parameter (test=1.e-5_R8)

    read (unit_struct,1000) title
    read (unit_struct,1010) lattic,nat,cform,irel
    REL=.TRUE.
    IF(IREL.EQ.'NREL') REL=.FALSE.
    ALLOCATE(aname(nat),mult(nat),jrj(nat),r0(nat),dx(nat),rmt(nat))
    allocate(zz(nat),rotloc(3,3,nat),v(nat))
    v=0.0d0
    ALLOCATE (pos(3,48*nat))
    read (unit_struct,1020) aa,bb,cc,alpha(1),alpha(2),alpha(3)
    IF(ABS(ALPHA(1)).LT.test) ALPHA(1)=90
    IF(ABS(ALPHA(2)).LT.test) ALPHA(2)=90
    IF(ABS(ALPHA(3)).LT.test) ALPHA(3)=90
    INDEX=0
    DO jatom=1,NAT
     INDEX=INDEX+1
     read(unit_struct,1030,iostat=ios) &
          iatnr, (pos(j,index),j=1,3), mult(jatom), isplit
     IF(ios /= 0) THEN
      WRITE(unit_out,*) iatnr,(pos(j,index),j=1,3),mult(jatom)
      WRITE(unit_out,*) 'ERROR IN STRUCT FILE READ'
      STOP
      ENDIF
      IF (mult(jatom) .EQ. 0) THEN
       WRITE (unit_out,6000) jatom, index, mult(jatom)
       STOP
       ENDIF
       DO m=1,mult(jatom)-1
        index=index+1
        READ(unit_struct,1031)iatnr,(pos(j,index),j=1,3)
       ENDDO
       READ(unit_struct,1050) aname(jatom),jrj(jatom),r0(jatom),rmt(jatom), &
            zz(jatom)
       dx(jatom)=LOG(rmt(jatom)/r0(jatom)) / (jrj(jatom)-1)
       rmt(jatom)=r0(jatom)*EXP( dx(jatom)*(jrj(jatom)-1) )
       READ(unit_struct,1051) ((rotloc(i,j,jatom),i=1,3),j=1,3)
    ENDDO
    ndif=index
    call realloc(pos, (/ 3, ndif /))
    ALLOCATE(rotij(3,3,ndif),tauij(3,ndif))
    READ(unit_struct,1151)iord
    ALLOCATE(iz(3,3,iord),tau(3,iord))
    DO j=1,iord
       READ(unit_struct,1101)((iz(j1,j2,j),j1=1,3),tau(j2,j),j2=1,3),inum
    ENDDO

1000 FORMAT(A80)
1010 FORMAT(A4,23X,I3,1x,a4,/,13X,A4,18X,A4)
1020 FORMAT(6F10.7,10X,F10.7)
1030 FORMAT(4X,I4,4X,F10.7,3X,F10.7,3X,F10.7,/,15X,I2,17X,I2)
1031 FORMAT(4X,I4,4X,F10.7,3X,F10.7,3X,F10.7)
1050 FORMAT(A10,5X,I5,5X,F10.9,5X,F10.5,5X,F10.5)
1051 FORMAT(20X,3F10.8)
1101 FORMAT(3(3I2,F11.8/),I8)
1151 FORMAT(I4)
6000 FORMAT(///,3X,'ERROR IN LAPW0 : MULT(JATOM)=0 ...', &
          /, 20X,'JATOM=',I3,3X,'INDEX=',I3,3X,'MULT=',I3)
  END SUBROUTINE init_struct

END MODULE struct

MODULE bessel
  use const, only: R8
  implicit none
  private :: R8

  real(R8), allocatable  ::rj(:,:),ri_mat(:,:)

CONTAINS
  subroutine init_bessel(LMAX2,LJMAX,NRAD,NRF)
    integer, intent(in) :: Lmax2, LJmax, Nrad, NRF

    integer :: l, l1, l2, lj

    L=0
    DO L1=0,LMAX2
       DO L2=0,LMAX2
          DO LJ=0,LJMAX
             IF (MOD((L1+L2+LJ),2) .EQ. 1) cycle
             IF ((L1+L2-LJ) .LT. 0)        cycle
             IF ((L1-L2+LJ) .LT. 0)        cycle
             IF ((-L1+L2+LJ) .LT. 0)       cycle
             L=L+1
          END DO
       END DO
    END DO
    ALLOCATE(rj(0:LJMAX+1,NRAD),ri_mat(NRF*NRF,L))
  END SUBROUTINE init_bessel
END MODULE bessel

MODULE xa
  use const, only: R8, C16
  implicit none
  private :: R8, C16
  public

  complex(C16), allocatable :: phs(:)
  real(R8),     allocatable :: fj(:,:),dfj(:,:),r(:)
  real(R8)                  :: bk(3),bkrot(3),bkrloc(3)

 CONTAINS
   subroutine init_xa(LMAX2,NMAT,NRAD,NB)
     integer, intent(in) :: Lmax2, Nmat, Nrad, NB

     allocate(phs(nb),fj(0:lmax2,nmat),dfj(0:lmax2,nmat))
     allocate(r(nrad))
   end subroutine init_xa
END MODULE xa

MODULE xa3
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

CONTAINS
  subroutine init_xa3(nb, nmat, num_kpts)
    integer, intent(in) :: NB, Nmat, num_kpts

    allocate( a  (Nmat, Nb, num_kpts), vecsz(num_kpts) )
    allocate( BKX(Nmat),         BKY(Nmat),         BKZ(Nmat),          &
         &    GX (Nmat,num_kpts), GY(Nmat,num_kpts), GZ(Nmat,num_kpts), &
         &    XK (     num_kpts), YK(     num_kpts), ZK(     num_kpts))

    GX = 0; GY = 0; GZ = 0
  END SUBROUTINE init_xa3
END MODULE xa3

module Amn_Mmn
  use const, only: C16

  implicit none
  private
  public :: overlap, init_Amn_Mmn, c

  complex(C16), allocatable :: overlap(:,:,:), c(:,:,:)

contains
  subroutine init_Amn_Mmn(nbands,npair,nproj)
    use param, only: Lmax2, Ndif

    integer, intent(in) :: nbands,npair,nproj

    allocate(overlap(nbands,nbands,npair), c(nproj,(LMAX2+1)*(LMAX2+1),ndif))

    overlap = 0; c = 0
  end subroutine init_Amn_Mmn
end module Amn_Mmn

module pairs
  implicit none
  public

  integer, allocatable :: KP(:), KPB(:)
  integer, allocatable :: BQX(:),BQY(:),BQZ(:), BQX1(:),BQY1(:),BQZ1(:)

contains
  subroutine  init_pairs(Npair)
    integer, intent(in) :: Npair

    allocate( kp  (Npair), kpb (Npair), &
         &    bqx (Npair), bqy (Npair), bqz (Npair), &
         &    bqx1(Npair), bqy1(Npair), bqz1(Npair))
  end subroutine init_pairs
end module pairs

module ams
  use const, only: R8

  real(R8), parameter :: atom_mass(103) = (/ &
         1.0,   4.0,   6.9,   9.0,  10.8,  12.0,  14.0,  16.0,  19.0,  20.2,&
        23.0,  24.3,  27.0,  28.1,  31.0,  32.0,  35.4,  40.0,  39.1,  40.0,&
        45.0,  47.9,  50.9,  52.0,  54.9,  55.8,  58.9,  58.7,  63.5,  65.4,&
        69.7,  72.6,  74.9,  79.0,  79.9,  83.8,  85.5,  87.6,  88.9,  91.2,&
        92.9,  95.9,  98.0, 101.1, 102.9, 106.4, 107.9, 112.4, 114.8, 118.7,&
       121.8, 127.6, 126.9, 131.3, 132.9, 137.3, 138.9, 140.1, 140.9, 144.2,&
       145.0, 150.4, 152.0, 157.3, 158.9, 162.5, 164.9, 167.3, 168.9, 173.0,&
       175.0, 178.5, 180.9, 183.8, 186.2, 190.2, 192.2, 195.1, 197.0, 200.6,&
       204.4, 207.2, 209.0, 209.0, 210.0, 222.0, 223.0, 226.0, 227.0, 232.0,&
       231.0, 238.0, 237.0, 244.0, 243.0, 247.0, 247.0, 251.0, 252.0, 257.0,&
       258.0, 259.0, 262.0 /)
end module ams

module lolog
  use param, only: Lmax2, lomax

  implicit none
  private

  integer, public :: nlo, nlov, nlon, n_rad(0:lmax2), ilo(0:lomax)
  logical, public :: loor(0:lomax), lapw(0:lmax2)
end module lolog

module gener
  use const, only: R8

  implicit none
  private

  ! transformation between u.c. and cartesian coordinates
  real(R8), public :: BR1(3,3), BR2(3,3)
end module gener

module atspdt
  use param, only: Lmax2, Nrf
  use const, only: R8

  implicit none
  private

  ! radial function and its slope at RMT
  real(R8), public :: P(0:Lmax2, Nrf), DP(0:Lmax2,Nrf)
end module atspdt

module loabc
  use param, only: lomax, Nloat, Nrf
  use const, only: R8

  implicit none
  private

  ! abc calculates the cofficients a,b,c of the lo
  real(R8), public :: alo(0:lomax, Nloat, Nrf)
end module loabc

module uhelp
  use param, only: Nrad
  use const, only: R8

  implicit none
  private

  real(R8), public :: A(NRAD), B(NRAD)
end module uhelp

module radfu
  use param, only: Nrad, Lmax2, Nrf
  use const, only: R8

  implicit none
  private

  ! radial functions large and small component
  real(R8), public :: RF1(NRAD,0:LMAX2,nrf), RF2(NRAD,0:LMAX2,nrf)
end module radfu

module rinteg
  use const, only: R8
  use param, only: Lmax2, Nrf

  implicit none
  private

  real(R8), public :: ri_mat(0:Lmax2, Nrf, Nrf)
end module rinteg

module PS1
  use const, only: R8
  implicit none
  private

  ! DEP,DEQ DERIVEES DE DP ET DQ   DB=ENERGIE/DVC    DVC VITESSE DE LA
  ! LUMIERE EN U.A.   DSAL=2.*DVC   DK NOMBRE QUANTIQUE KAPPA
  ! DM=PAS EXPONENTIEL/720., DKOEF=1./720.
  real(R8), public :: DEP(5), DEQ(5), DB, DVC, DSAL, DK, DM
end module PS1

!!/---
!! Local Variables:
!! mode: f90
!! End:
!!\---
!!

!! Time-stamp: <2016-07-06 12:29:24 assman@faepop71.tu-graz.ac.at>
