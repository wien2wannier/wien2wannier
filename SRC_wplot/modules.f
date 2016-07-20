!!! wien2wannier/SRC_wplot/modules.f

module wplot
  use const, only: BUFSZ

  implicit none

  private :: BUFSZ
  public

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

  integer, parameter :: unit_grid=7,  unit_psink=21, unit_psiarg=22
  integer, parameter :: unit_inwf=31, unit_chk  =32, unit_rot   =33

  integer, parameter :: Nsym=48, Lmax7=8

!!! The following is set by ‘wplot.f’ to be passed to ‘main.F’
  character(BUFSZ) :: vecfn, psinkfn, psiargfn, outfn
  integer          :: idx_wann=0, iproc=0
end module wplot


!---------------- Reading ‘chk’ files                  ----------------------
!
! Mostly copied from Wannier90
module Wannier90
  use const, only: DPk
  use util,  only: newunit

  implicit none
  private
  public :: chk_t, chk_read

  interface chk_read
     module procedure chk_read_fname, chk_read_unit, chk_read_argstr
  end interface chk_read

  type chk_t
     character(len=33)          :: header ! usually, date and time
     character(len=20)          :: checkpoint
     !                             #k,       #WF       #B
     integer                    :: num_kpts, num_wann, num_bands
     ! total number of neighbors per k-point
     integer                    :: nntot

     ! selection of the optimal subspace (#B × #WF × #k)
     complex(DPk), allocatable  :: u_matrix_opt(:,:,:)
     ! rotation to the optimally smooth states in subspace (#WF × #WF × #k)
     complex(DPk), pointer      :: u_matrix    (:,:,:)
     ! overlaps M_mn are stored in ‘chk’ for restart (#WF × #WF × #nn × #k)
     complex(DPk), allocatable  :: m_matrix(:,:,:,:)
     !                            [Å] (3 × #WF),        [Å²] (#WF)
     real(DPk),    allocatable :: wannier_centres(:,:), wannier_spreads(:)
     ! kpoints in lattice vecs (3 × #k)
     real(DPk),    allocatable :: kpt_latt(:,:)

     ! disentanglement parameters
     logical                    :: have_disentangled
     real(DPk)                  :: omega_invariant
     logical, allocatable       :: lwindow(:,:) ! (#B × #k)
     integer, allocatable       :: ndimwin(:)   ! (#k)

     integer                    :: num_exclude_bands
     integer, allocatable       :: exclude_bands(:) ! (num_exclude_bands)

     real(DPk), dimension(3, 3) :: real_lattice, recip_lattice
     integer,   dimension(3)    :: mp_grid
  end type chk_t

contains
  subroutine chk_read_unit(lun, chk, read_mmn)
    integer,           intent(in)  :: lun
    type(chk_t),       intent(out) :: chk
    logical, optional, intent(in)  :: read_mmn

    logical   :: rm
    integer   :: j, jk, l

    if (present(read_mmn)) then
       rm = read_mmn
    else
       rm = .true.
    end if

    read(lun) chk%header
    read(lun) chk%num_bands
    read(lun) chk%num_exclude_bands
    allocate( chk%exclude_bands(chk%num_exclude_bands) )
    read(lun) chk%exclude_bands
    read(lun)(chk%real_lattice (:,j), j=1,3)
    read(lun)(chk%recip_lattice(:,j), j=1,3)
    read(lun) chk%num_kpts
    read(lun) chk%mp_grid
    allocate( chk%kpt_latt(3, chk%num_kpts) )
    read(lun)(chk%kpt_latt(:,j), j=1,chk%num_kpts)
    read(lun) chk%nntot
    read(lun) chk%num_wann
    read(lun) chk%checkpoint
    read(lun) chk%have_disentangled

    DIS: if (chk%have_disentangled) then
       read(lun)  chk%omega_invariant
       allocate(  chk%lwindow(chk%num_bands, chk%num_kpts) )
       read(lun) (chk%lwindow(:,jk), jk=1,chk%num_kpts)
       allocate(  chk%ndimwin(chk%num_kpts) )
       read(lun)  chk%ndimwin
       allocate(  chk%u_matrix_opt(chk%num_bands, chk%num_wann, chk%num_kpts) )
       read(lun)((chk%u_matrix_opt(:,j,jk),j=1,chk%num_wann),jk=1,chk%num_kpts)
    end if DIS

    allocate(  chk%u_matrix(chk%num_wann, chk%num_wann, chk%num_kpts) )
    read(lun)((chk%u_matrix(:,j,jk), j=1,chk%num_wann), jk=1,chk%num_kpts)

    Mmn: if (rm) then
       allocate(chk%m_matrix(chk%num_wann,chk%num_wann,chk%nntot,chk%num_kpts))
       read(lun) (((chk%m_matrix(:, j, l, jk), &
            &            j=1,chk%num_wann), l=1,chk%nntot), jk=1,chk%num_kpts)
    else
       read(lun)
    end if Mmn

    allocate( chk%wannier_centres(3, chk%num_wann) )
    read(lun)(chk%wannier_centres(:,j),j=1,chk%num_wann)
    allocate( chk%wannier_spreads(chk%num_wann) )
    read(lun) chk%wannier_spreads
  end subroutine chk_read_unit


!!! Wrappers for calling the above with a filename instead of an open
!!! unit
  subroutine chk_read_fname(fname, chk, read_mmn)
    use util, only: newunit

    character(*),      intent(in)  :: fname
    type(chk_t),       intent(out) :: chk
    logical, optional, intent(in)  :: read_mmn

    integer :: lun
    logical :: rm

    if (present(read_mmn)) then
       rm = read_mmn
    else
       rm = .true.
    end if

    open(unit=newunit(lun), file=fname, status='old', FORM='unformatted')
    call chk_read_unit(lun, chk, rm)
    close(lun)
  end subroutine chk_read_fname

  subroutine chk_read_argstr(arg, chk, read_mmn)
    use clio, only: argstr

    type(argstr),      intent(in)  :: arg
    type(chk_t),       intent(out) :: chk
    logical, optional, intent(in)  :: read_mmn

    logical :: rm

    if (present(read_mmn)) then
       rm = read_mmn
    else
       rm = .true.
    end if

    call chk_read_fname(arg%s, chk, rm)
  end subroutine chk_read_argstr
end module Wannier90


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

module PS1
  use const, only: R8
  implicit none
  private

  ! DEP,DEQ DERIVEES DE DP ET DQ   DB=ENERGIE/DVC    DVC VITESSE DE LA
  ! LUMIERE EN U.A.   DSAL=2.*DVC   DK NOMBRE QUANTIQUE KAPPA
  ! DM=PAS EXPONENTIEL/720., DKOEF=1./720.
  real(R8), public :: DEP(5), DEQ(5), DB, DVC, DSAL, DK, DM
end module PS1

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

  real(DPk) :: TRANS(3,NSYM)
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

module uhelp
  use param, only: Nrad
  use const, only: R8

  implicit none
  private

  real(R8), public :: A(NRAD), B(NRAD)
end module uhelp

!!/---
!! Local Variables:
!! mode: f90
!! End:
!!\---
!!
!! Time-stamp: <2016-07-20 11:38:09 assman@faepop71.tu-graz.ac.at>
