!!! wien2wannier/SRC_wplot/modules.f
!!!
!!!    Modules for wplot.  This file contains modules that are
!!!    independent of real/complex compilation.

!!/=== Defined here: =============================
!!
!! wplot:     AddLoc, kconjg, mvatom, UseRot, Nsym, Lmax7, idx_wann, iproc,
!!            unit_chk, unit_inwf, unit_grid, unit_psiarg, unit_psink,
!!            unit_rot, vecfn, psinkfn, psiargfn, outfn
!!
!! Wannier90: type(chk_t), chk_read()
!!
!! atspdt:    P(:,:,:), DP(:,:,:)
!!
!! bessfu:    FJ(:,:,:), DFJ(:,:,:), IRAD(:), RAD(:)
!!
!! grid:      ilat(:,:), ireg(:), iri(:)
!!
!! latt:      BR1(3,3), BR2(3,3), BR3(3,3), BR4(3,3)
!!
!! loabc:     ALO(:,:,:,:)
!!
!! lolog:     ILO(:,:), LAPW(:,:), NLO
!!
!! radfu:     RRAD(:,:,:,:)
!!
!! struct:    POS(:,:), RMT(:)
!!
!! sym2:      imat(3,3,Nsym), iord, trans(3,Nsym)
!!
!! work:      aug(:,:,:)
!!
!! work1:     A(Nrad), B(Nrad)
!!
!! Procedure modules:
!!
!!    auggen_m, auglo_m, augpw_m, bessel_m, findmt_m, gbass_m,
!!    latgen_m, locdef_m, rint13_m, rotdef_m, spcgen_m, trans_m,
!!    wavint_m, wavsph_m
!!
!!\===============================================

module wplot
  use const, only: BUFSZ

  implicit none

  private :: BUFSZ
  public; save

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
  character(BUFSZ) :: vecfn, psinkfn, psiargfn, outfn, gridfn
  integer          :: idx_wann=0, iproc=0
end module wplot


module latt
  use const, only: DPk

  implicit none
  private; save

! BR1(i,:) -- the real space lattice vectors a_i of the conventional u.c.
! BR2(i,:) -- the real space lattice vectors a_i of the primitive unit cell
! BR3(i,:) -- the reciprocal lattice vectors b_i of the conventional u.c
! BR4(i,:) -- the reciprocal lattice vectors b_i of the primitive unit cell

  real(DPk), public :: BR1(3,3), BR2(3,3), BR3(3,3), BR4(3,3)
end module latt

!---------------- Reading ‘inwplot’ files              ----------------------
module inwplotmod
  use const, only: DPk

  implicit none
  private
  public :: inwplot_t, inwplot_read, unit_ATU, unit_ANG

  interface inwplot_read
     module procedure &
          inwplot_read_fname, inwplot_read_unit, inwplot_read_argstr
  end interface inwplot_read

  real(DPk), parameter :: unit_ATU=1, unit_ANG=1/sqrt(0.529177_DPk**3)

  type inwplot_t
     logical      :: checkortho, dephas, WFrot, large, extgrid
     integer      :: WFidx
     real(DPk)    :: unit
     character(3) :: unit_name
     real(DPk)    :: orig(3), xend(3), yend(3), zend(3)

     real(DPk), allocatable :: rgrid(:,:)
  end type inwplot_t

contains

  subroutine inwplot_read_unit(lun, inwplot)
!!! This is based on grdgen.f (including some comments)
    use const,  only: ORTHO_TEST, TAU
    use param,  only: unit_out
    use util,   only: uppercase, string
    use clio,   only: croak
    use wplot,  only: unit_grid, unit_psink, gridfn
    use latt,   only: br1, br4

    integer,         intent(in)  :: lun
    type(inwplot_t), intent(out) :: inwplot

    character(len=5) :: relcomp
    character(len=3) :: mode, post, unit
    character(len=1) :: flag
    integer          :: Npt, ig, div, Npx, Npy, Npz, i, j, ix,iy,iz
    real(DPk)        :: cartax(3,3), fracax(3,3), prodax(3,3), angles(3)
    real(DPk)        :: cartorig(3), fracorig(3), cosphi, phi(3,3)
    character        :: xyz(3) = (/ 'x', 'y', 'z' /)

    character(*), parameter :: &
         fmt_orig  = "(' origin', 1X, 3F9.5, 3X, 3F13.7)",   &
         fmt_axend = "(' ',A1,'-axis',1X,3F9.5,3X,3F13.7)",  &
         fmt_axabs = "(' ‖',A1,'‖ ',F13.7)",                 &
         fmt_angle = "(' ∠(',A1,',',A1,') ',F10.7,3X,F8.3)", &
         fmt_grid  = "('#', i7, 2x, f13.7, 2(2x, f10.5))"

    read(lun, '(A3, A1)') mode, flag

    gridgen: select case ( uppercase(adjustl(mode)) )
    case ('ANY')
       inwplot%extgrid = .true.
       read(lun, *) Npt
       allocate(inwplot%rgrid(3, Npt))

       cartfrac: select case (uppercase(flag))
       case ('C')
          do ig=1,Npt
             read(unit_grid, *) inwplot%rgrid(:, ig), div
             inwplot%rgrid(:, ig) = inwplot%rgrid(:, ig) / div
          end do

       case ('F')
          do ig=1,Npt
             read(unit_grid, *) inwplot%rgrid(:, ig)
             inwplot%rgrid(:, ig) = matmul(br4, inwplot%rgrid(:, ig))
          end do

       case default
          call croak("unknown FLAG `"//flag &
            &   //"' in inwplot_read() [should be `C' or `F']")
       end select cartfrac

       write(unit_psink,                                  &
            "('ANY   NP =',I0,/,                          &
            & 'order according to the grid points ',/,    &
            & 'provided in the input file ',/,            &
            & '`', A, '''')") &
            Npt, trim(gridfn)

       write(unit_out, &
            & "(' ARBITRARY LIST OF GRID POINTS' / &
            & ' -----------------------------' / &
            & ' number of grid points: ',I7)") Npt

    case ('3D')
       inwplot%extgrid = .false.
       select case (uppercase(flag))
       case ('O', ' ')
          inwplot%checkortho = .true.
       case ('N')
          inwplot%checkortho = .false.
       case default
          call croak("unknown FLAG `"//flag &
          &        //"' in inwplot_read() [should be `O' or `N']")
       end select

       read(lun, *) inwplot%orig, div; inwplot%orig = inwplot%orig/div
       read(lun, *) inwplot%xend, div; inwplot%xend = inwplot%xend/div
       read(lun, *) inwplot%yend, div; inwplot%yend = inwplot%yend/div
       read(lun, *) inwplot%zend, div; inwplot%zend = inwplot%zend/div

       read(lun, *) Npx, Npy, Npz

       allocate (inwplot%rgrid(3, Npx*Npy*Npz))

       write(unit_out, &
            &"(' 3D-NET OF GRID POINTS'                                    /&
            &' ---------------------'                                      /&
            &' number of grid points for x, y, z: ',3(I0,' '),              &
            &' (total: ',I0,')')")                                          &
            Npx, Npy, Npz, Npx*Npy*Npz

       write(unit_out, "(/' PLOTTING AREA'                         &
            &            /' -------------'                         &
            &            /' x = Sum(j=1,3) f_i a_i  with  f_i in', &
            &             ' conventional fractional coordinates'   &
            &           //'           f_1      f_2      f_3   ',   &
            &             '        x [Bohr]     y [Bohr]     z [Bohr]')")

       !---------------------------------------------------------------------
       ! Transform origin and axes into primitive fractional coordinates
       ! a) transform into Cartesian coordinates
       !    and check axes for orthogonality
       fracax(:,1) = inwplot%xend-inwplot%orig
       fracax(:,2) = inwplot%yend-inwplot%orig
       fracax(:,3) = inwplot%zend-inwplot%orig

       cartax   = matmul(transpose(br1), fracax)
       cartorig = matmul(transpose(br1), inwplot%orig)
       prodax   = matmul(transpose(cartax), cartax)

       write(unit_out, fmt_orig) inwplot%orig, cartorig
       do i=1,3
          write(unit_out, fmt_axend) xyz(i), fracax(:,i), cartax(:,i)
       end do

       write(unit_out, "(/'         length [Bohr]')")
       do i=1,3
          write(unit_out, fmt_axabs) xyz(i), sqrt(prodax(i,i))
       end do

       write(unit_out, "(/'         cos(φ)        φ [°]')")
       phi = 0
       do i=1, 2
          do j=i+1, 3
             cosphi = prodax(i,j) / sqrt(prodax(i,i) * prodax(j,j))
             phi(i,j) = acos(cosphi) * 360 / TAU
             phi(j,i) = phi(i,j)

             write(unit_out, fmt_angle) xyz(i), xyz(j), cosphi, phi(i,j)
          end do
       end do

       angles = (/ phi(2,3), phi(1,3), phi(1,2) /)
       if (inwplot%checkortho &
            .and. any(abs(angles - 90) > ORTHO_TEST/360*TAU)) &
            then
          call croak ('nonorthogonal axes (angles [°]: ' &
               &      // trim(string(angles)) // ')')
       end if

       ! b) transform into primitive fractional coordinates
       fracax   = matmul(br4, cartax)
       fracorig = matmul(br4, cartorig)
       !---------------------------------------------------------------------

       ! Generate evaluation grid
       ! use gnuplot order for data generation, i.e.
       ! POS = IZ + (IY-1)*NP(3) + (IX-1)*NP(3)*NP(2)
       fracax(:,1) = fracax(:,1) / max(Npx-1, 1)
       fracax(:,2) = fracax(:,2) / max(Npy-1, 1)
       fracax(:,3) = fracax(:,3) / max(Npz-1, 1)

       forall (ix=0:Npx-1, iy=0:Npy-1, iz=0:Npz-1)
          inwplot%rgrid(:, 1 + iz + + Npz*iy + Npz*Npy*ix) = &
               fracorig + ix*fracax(:,1) + iy*fracax(:,2) + iz*fracax(:,3)
       end forall

       write(0, *) inwplot%rgrid

       write(unit_psink, "('# 3D  NP      length(α)    ∠(α, x)    ∠(α, y)')")

       write(unit_psink, fmt_grid) Npx, sqrt(prodax(1,1))
       write(unit_psink, fmt_grid) Npy, sqrt(prodax(2,2)), phi(1,2)
       write(unit_psink, fmt_grid) Npz, sqrt(prodax(3,3)), phi(1,3), phi(2,3)

       write(unit_psink, &
            "('#order: (((w(x,y,z), z=1,#z), y=1,#y), x=1,#x)')")
       stop

    case default
       call croak("unknown MODE `"//mode &
            &   //"' in inwplot_read() [should be `ANY' or `3D']")
    end select gridgen

    read(lun, '(A3)') post
    select case ( uppercase(adjustl(post)) )
    case ('DEP'); inwplot%dephas = .true.
    case ('NO');  inwplot%dephas = .false.
    case default; call croak("unknown postprocessing option `"//post//"'")
    end select

    read(lun, *) unit, relcomp
    select case ( uppercase(adjustl(unit)) )
    case ('', 'AU', 'ATU');
       inwplot%unit = unit_ATU
       inwplot%unit_name ='ATU'
    case ('ANG');
       inwplot%unit = unit_ANG
       inwplot%unit_name ='ANG'
    case default
       call croak("unknown units option `"//unit//"'")
    end select
    select case ( uppercase(adjustl(relcomp)) )
    case ('LARGE', ''); inwplot%large = .true.
    case ('SMALL');     inwplot%large = .false.
    case default
       call croak("unknown relativistic component `"//relcomp//"'")
    end select

    read(lun, *) inwplot%WFidx, inwplot%WFrot
  end subroutine inwplot_read_unit

  subroutine inwplot_read_fname(fname, inwplot)
    use util, only: newunit

    character(*),    intent(in)  :: fname
    type(inwplot_t), intent(out) :: inwplot

    integer :: lun

    open(unit=newunit(lun), file=fname, status='old')
    call inwplot_read_unit(lun, inwplot)
    close(lun)
  end subroutine inwplot_read_fname

  subroutine inwplot_read_argstr(arg, inwplot)
    use clio, only: argstr

    type(argstr),    intent(in)  :: arg
    type(inwplot_t), intent(out) :: inwplot

    call inwplot_read_fname(arg%s, inwplot)
  end subroutine inwplot_read_argstr
end module inwplotmod


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
  private; save

  real(DPk), allocatable, public :: P(:,:,:), DP(:,:,:)
end module atspdt

module bessfu
  use const, only: DPk

  implicit none
  private; save

  real(DPk), allocatable, public :: FJ(:,:,:), DFJ(:,:,:), RAD(:)
  integer,   allocatable, public :: IRAD(:)
end module bessfu

module grid
  use const, only: DPk

  implicit none
  private; save

  integer,   allocatable, public :: ireg(:),ilat(:,:),iri(:)
end module grid

module loabc
  use const, only: DPk

  implicit none
  private; save

  real(DPk), allocatable, public :: ALO(:,:,:,:)
end module loabc

module lolog
  implicit none
  public; save

  integer              :: NLO
  logical, allocatable :: LAPW(:,:)
  integer, allocatable :: ILO(:,:)
end module lolog

module radfu
  use const, only: DPk

  implicit none
  private; save

  real(DPk), allocatable, public :: RRAD(:,:,:,:)
end module radfu

module sym2
  use const, only: DPk
  use wplot, only: Nsym

  implicit none
  private; save

  real(DPk), public :: TRANS(3,NSYM)
  integer,   public :: IMAT(3,3,NSYM), IORD
end module sym2

module work
  use const, only: DPK

  implicit none
  private; save

  complex(DPK), allocatable, public :: aug(:,:,:)
end module work


!---------------------------  Procedure modules  ---------------------------
module     rint13_m; contains
subroutine rint13(rel, A, B, X, Y, S, jatom, stru)
  use const,     only: DPk, clight
  use param,     only: Nrad
  use structmod, only: struct_t

  implicit none

  logical,        intent(in)  :: rel
  real(DPk),      intent(in)  :: A(Nrad), B(Nrad), X(Nrad), Y(Nrad)
  real(DPk),      intent(out) :: S
  integer,        intent(in)  :: jatom
  type(struct_t), intent(in)  :: stru

  real(DPk) :: cin, d, R, R1, Z2, Z4, P1, P2
  integer   :: j, j1

!     last changes: 01.11.00 ub (updating comments)
!
!     PERFORM RADIAL INTEGRALS
!     Int(0,Rmt) u(r) v(r) + 1/c^2 us(r) vs(r) r^2 dr
!
!     with c = 137.037 (274.074) in Hartree (Rydberg) units
!
!     For non-relativistic calculations c := 10^+11 Rydberg units.
!
!----------------------------------------------------------------------------
! Input:
! REL    .TRUE. for scalar relativistc calculations
! A(:)   r * u (r) on the radial mesh
! B(:)   r * us(r) on the radial mesh
! X(:)   r * v (r) on the radial mesh
! Y(:)   r * vs(r) on the radial mesh
! JATOM  the current type of atom
!
! from RADGRD
! DX  (JATOM)  logaritmic increment of the radial mesh
!
! Output:
! S    the value of the radial integral
!----------------------------------------------------------------------------

  CIN = 1/clight**2
  if(.not.REL) CIN=4E-22    ! legacy value

!     CONVERT FROM RYDBERG TO HARTREE UNITS
  D=exp(stru%dx(JATOM))
  J=3-mod(STRU%NPT(JATOM),2)
  J1=J-1
  R=STRU%R0(JATOM)*(D**(J-1))
  R1=R/D
  Z4=0
  Z2=0
10 Z4=Z4+R*(A(J)*X(J)+CIN*B(J)*Y(J))
  R=R*D
  J=J+1
  if(J.ge.STRU%NPT(JATOM)) goto 20
  Z2=Z2+R*(A(J)*X(J)+CIN*B(J)*Y(J))
  R=R*D
  J=J+1
  goto 10
20 P1=STRU%R0(JATOM)*(A(1)*X(1)+CIN*B(1)*Y(1))
  P2=R1*(A(J1)*X(J1)+CIN*B(J1)*Y(J1))
  S=2*Z2+4*Z4+R*(A(J)*X(J)+CIN*B(J)*Y(J))+P2
  S=(stru%dx(JATOM)*S+P1)/3.0D0
  if(J1.gt.1) S=S+0.5D0*stru%dx(JATOM)*(P1+P2)
end subroutine rint13
end module     rint13_m


module     gbass_m; contains
subroutine gbass(RBAS,GBAS,TWOPI)
  use const, only: DPk, TAU

  implicit none

  real(DPk), dimension(3,3) :: RBAS, GBAS
  logical                   :: TWOPI

  intent(in)  :: Rbas, twoPI
  intent(out) :: Gbas

  ! << Input >>
  ! RBAS(i,:) -- the real space lattice vectros a_i
  ! TWOPI     -- normalization of the reciprocal lattice vectors:
  !                if .TRUE.  let < a_i | b_j > = 2pi * delta(i,j)
  !                if .FALSE. let < a_i | b_j > = delta(i,j)
  !
  ! << Output >>
  ! GBAS(j,:) -- the corresponding reciprocal lattice vectors b_j

  real(DPk) :: Vuc, fac
  integer   :: i

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
  VUC=0
  do I=1,3
     VUC = VUC + RBAS(1,I)*GBAS(1,I)
  end do

  if(TWOPI)then
     FAC = TAU/2 / VUC
  else
     FAC = 1  / VUC
  endif

  GBAS = GBAS * FAC
end subroutine gbass
end module     gbass_m


module     spcgen_m; contains
subroutine spcgen(NAT,RMT,ATMS)
  use const, only: DPk
  use latt,  only: br4

  implicit none

  integer,   intent(in)  :: NAt
  real(DPk), intent(in)  :: RMT(NAT)
  real(DPk), intent(out) :: ATMS(3,NAT)

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

  real(DPk) :: B(3)
  integer   :: i, jatom

  do I=1,3
     B(I) = sqrt( BR4(I,1)**2 + BR4(I,2)**2 + BR4(I,3)**2 )
  end do
  do JATOM=1,NAT
     do I=1,3
        ATMS(I,JATOM) = RMT(JATOM) * B(I)
     end do
  end do
end subroutine spcgen
end module     spcgen_m


module     wavint_m; contains
subroutine wavint(R,NPW,PSI,bk,coef,nmat)
  use const, only: DPk

  implicit none

  integer,      intent(in)  :: NPW, NMat
  real(DPk),    intent(in)  :: R(3), BK(3,Nmat)
  complex(DPk), intent(out) :: Psi
  complex(DPk), intent(in)  :: coef(NMat)

! evaluation of the wave function in the interstitial region:
! psi(r) = Sum(K) c_K/sqrt(V) e^i(K+k)r
! --------------------------------------------------------------
! Input:
! R    -- grid point in (gobal) Cartesian coordinates
! NPW  -- number of PW basis functions
!
! MODULE EIGVEC
! BK   -- the PW wave vectors K+k in (gobal) Cartesian coordinates
! COEF -- the PW coefficients c_K (including 1/sqrt(V))
!
! Output:
! PSI  -- the wavr function psi(r)
! --------------------------------------------------------------

  real(DPk) :: arg
  integer   :: iPW

  PSI = 0
  do IPW=1,NPW
     ARG = R(1)*BK(1,IPW) + R(2)*BK(2,IPW) + R(3)*BK(3,IPW)
     PSI = PSI + COEF(IPW) * cmplx(cos(arg), sin(arg), DPk)
  end do
end subroutine wavint
end module     wavint_m


module     auggen_m; contains
subroutine auggen(stru, large)
  use lolog,     only: Nlo, ilo, lapw
  use loabc,     only: Alo
  use atspdt,    only: P, DP
  use radfu,     only: rrad
  use uhelp,     only: A, B
  use const,     only: DPk
  use param,     only: NLOat, Nrad, Nrf, LOmax, &
       &               unit_out, unit_vsp, unit_vector
  use wplot,     only: Lmax7
  use structmod, only: struct_t, radgrid

  !! procedure includes
  use diracout_m
  use outwin_m
  use rint13_m
  use dergl_m

  implicit none

  type(struct_t), intent(in) :: stru
  logical,        intent(in) :: large

  logical   :: rlo(1:nloat,0:lomax)
  integer   :: i, k, l, m, jatom, jlo, jrf, irf, nodes, kappa
  real(DPk) :: cfac, RMT2, rnorm, uve, duve, uv, uvb, duvb, duv
  real(DPk) :: dele, delei, fl, ei, e1, cross, clight, r_m
  real(DPk) :: pi12lo, pe12lo, xac, xbc, xcc, alonorm
  real(DPk) :: E(0:LMAX7), ELO(0:LOMAX,NLOAT), PEI(0:LMAX7)

  real(DPk), dimension(Nrad)             :: AE, BE, VR
  real(DPk), dimension(NRAD,0:LMAX7,NRF) :: RAD1, RAD2

  cfac = 1 / 274.074_DPk
  if(.not. stru%rel) cfac = 1e-11_DPk
  if(large) cfac = 1

  ! << skip header of *.vsp file >>
  read(unit_vsp,1000)

  write(unit_out,2000)
  NLO = 0
  ALO = 0
  ILO = 0
  do JATOM=1,STRU%NNEQ
     RMT2 = stru%rmt(jatom)**2
     write(unit_out,2010) jatom, stru%RMT(jatom)

     ! << read atomic potential r * V(r) and convert into Rydberg >>
     read(unit_vsp,1010)
     read(unit_vsp,1020) VR(1 : stru%Npt(jatom))
     read(unit_vsp,1030)

     VR(1 : stru%Npt(jatom)) = VR(1 : stru%Npt(jatom)) / 2

     ! << read augmentation energies and check for local orbitals >>
     read(unit_vector) E
     read(unit_vector) ELO
     do L=0,LMAX7
        LAPW(L,JATOM)=.true.
        if (E(L) > 150) then
           E(L)=E(L)-200.d0
           LAPW(L,JATOM)=.false.
        end if
     end do

     do L=0,LOMAX
        do k=1,nloat
           rlo(k,l)=.false.
           if(ELO(L,K).lt.995D0)then
              ilo(L,jatom)=ilo(L,jatom)+1
              if(.not.lapw(l,jatom).and.k.eq.1) goto 666
              if(k.eq.nloat) then
                 rlo(ilo(l,jatom),l)=.true.
                 goto 666
              endif
666           continue
              NLO=NLO+(2*L+1)*stru%MULT(JATOM)
           end if
        end do
     end do

     !       << regular augmentation functions >>
     !       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     write(unit_out,2020)
     DELE=2.0D-3
     DELEI=0.25D0/DELE
     do l=0,LMAX7
        FL=L
        EI=E(L)/2.0d0
        !         << compute d/dE u_l(r,E) by finite differences >>
        E1=EI-DELE
        call outwin(stru, jatom, Vr, E1, FL, UVB, DUVB, Nodes)
        call RINT13(stru%rel,A,B,A,B,RNORM,JATOM, stru)
        RNORM = 1.0D0/sqrt(RNORM)
        do M=1,stru%Npt(jatom)
           AE(M) = RNORM * A(M)
           BE(M) = RNORM * B(M)
        end do
        UVB  = RNORM * UVB
        DUVB = RNORM * DUVB
        E1=EI+DELE
        call outwin(stru, jatom, Vr, E1, FL, UVE, DUVE, Nodes)
        call RINT13(stru%rel,A,B,A,B,RNORM,JATOM, stru)
        RNORM = 1.0D0/sqrt(RNORM)
        UVE  = DELEI*(RNORM*UVE -UVB )
        DUVE = DELEI*(RNORM*DUVE-DUVB)
        do M=1,stru%Npt(jatom)
           AE(M) = DELEI*(RNORM*A(M)-AE(M))
           BE(M) = DELEI*(RNORM*B(M)-BE(M))
        end do

        !         << now compute u_l(r,E) >>
        call outwin(stru, jatom, Vr, EI, FL, UV, DUV, Nodes)
        call RINT13(stru%rel,A,B,A,B,RNORM,JATOM, stru)
        RNORM = 1.0D0/sqrt(RNORM)
        do M=1,stru%Npt(jatom)
           A(M) = RNORM*A(M)
           B(M) = RNORM*B(M)
        end do
        P(L,1,JATOM) = RNORM*UV
        DP(L,1,JATOM) = RNORM*DUV
        !
        !         << insure orthogonality of d/dE u_l(r,E) on u_l(r,E) >>
        !
        call RINT13(stru%rel,A,B,AE,BE,CROSS,JATOM, stru)
        do M=1,stru%Npt(jatom)
           AE(M) = (AE(M)-CROSS*A(M))
           BE(M) = (BE(M)-CROSS*B(M))
        end do
        P(L,2,JATOM) = UVE -CROSS*P(L,1,JATOM)
        DP(L,2,JATOM) = DUVE-CROSS*DP(L,1,JATOM)
        do I=1,stru%Npt(jatom)
           RAD1(I,L,1) = A(I)
           RAD1(I,L,2) = AE(I)
           RAD2(I,L,1) = B(I)
           RAD2(I,L,2) = BE(I)
        end do
        call RINT13(stru%rel,AE,BE,AE,BE,PEI(L),JATOM, stru)
        write(unit_out,2030) L,E(L),P(L,1,JATOM),DP(L,1,JATOM),P(L,2,JATOM),DP(L,2,JATOM)
     end do
     !
     !       << local orbitals >>
     do L=0,LOMAX
        irf=2
        do jlo=1,ilo(L,jatom)
           if (lapw(l,jatom).or.(jlo.gt.1)) then
              irf=irf+1
              DELE=2.0D-3
              DELEI=0.25D0/DELE
              FL=L
              EI=elo(l,jlo)/2.d0
              if(rlo(jlo,l)) then
                 ei=elo(l,nloat)/2.d0
                 kappa=l
                 call diracout(stru, jatom, Vr, ei, kappa, uv, duv, nodes)
                 call dergl(stru, jatom, a, b)
                 do m=1,stru%Npt(jatom)
                    r_m=stru%R0(JATOM)*exp(stru%dx(jatom)*(m-1))
                    b(m)=b(m)*r_m/(2.d0*clight+(elo(l,jlo)- &
                         2.d0*vr(m)/r_m)/(2.d0*clight))
                    b(m)=b(m)*clight
                 enddo
              else
                 call outwin(stru, jatom, Vr, EI, FL, UV, DUV, Nodes)
              endif
              call RINT13(stru%rel,A,B,A,B,RNORM,JATOM, stru)
              RNORM = 1.0d0/sqrt(RNORM)
              do M=1,stru%Npt(jatom)
                 RAD1(M,L,irf) = RNORM*A(M)
                 RAD2(M,L,irf) = RNORM*B(M)
              enddo
              P(L,irf,jatom)  = RNORM*UV
              DP(L,irf,jatom) = RNORM*DUV
              call RINT13(stru%rel,RAD1(1,L,1),RAD2(1,L,1),RAD1(1,L,irf),RAD2(1,L,irf),PI12LO,JATOM, stru)
              call RINT13(stru%rel,RAD1(1,L,2),RAD2(1,L,2),RAD1(1,L,irf),RAD2(1,L,irf),PE12LO,JATOM, stru)
           endif

           if (LAPW(L,JATOM)) then
              XAC=(P(L,irf,jatom)*DP(L,2,JATOM)-DP(L,irf,jatom)*P (L,2,JATOM))*RMT2
              XBC=(P(L,1,jatom)*DP(L,irf,JATOM)-P(L,irf,JATOM)*DP(L,1,jatom))*RMT2
              XCC=XAC*(XAC+2.0D0*PI12LO) &
                   + XBC*(XBC*PEI(L)+2.0D0*PE12LO)+1.0D0
              ALO(L,jlo,irf,JATOM) = 1.0D0/max( sqrt(XCC) , 0.005D0 )
              ALO(L,jlo,1,JATOM)=XAC*ALO(L,jlo,irf,JATOM)
              ALO(L,jlo,2,JATOM)=XBC*ALO(L,jlo,irf,JATOM)
           else
              if (jlo.eq.1) then
                 alonorm=sqrt(1.d0+(P(L,1,JATOM)/P(L,2,JATOM))**2*PEI(L))
                 ALO(L,jlo,1,JATOM)=1.d0/alonorm
                 ALO(L,jlo,2,JATOM)=-P(L,1,JATOM)/P(L,2,JATOM)/alonorm
              else
                 xbc=-P(l,1,jatom)/P(L,irf,jatom)
                 xac=sqrt(1+xbc**2+2*xbc*PI12LO)
                 ALO(l,jlo,1,jatom)=1.d0/xac
                 ALO(l,jlo,irf,jatom)=xbc/xac
              endif
           endif
        end do
     end do
     !
     do L=0,LOMAX
        do JLO=1,ilo(l,jatom)
           write(unit_out,2070)L,(alo(l,jlo,jrf,jatom),jrf=1,nrf)
        end do
     end do

     !       << re-scale radial augmentation functions  >>
     !       << with VR just being a working array here >>
     Vr(1:stru%Npt(jatom)) = &
          cfac / radgrid(stru, jatom, (/( i, i=1,stru%Npt(jatom) )/))

     do L=0,LMAX7
        if(LARGE)then
           do irf=1,nrf
              do I=1,stru%Npt(jatom)
                 RRAD(I,L,irf,JATOM) = RAD1(I,L,irf) * VR(I)
              enddo
           enddo
        else
           do irf=1,nrf
              do I=1,stru%Npt(jatom)
                 RRAD(I,L,irf,JATOM) = RAD2(I,L,irf) * VR(I)
              enddo
           enddo
        endif
     enddo
  end do

!
 1000 format(//)
 1010 format(/////)
 1020 format(3X,4E19.12)
 1030 format(/////)
!
 2000 format(/' AUGMENTATION FUNCTIONS' &
             /' ----------------------' &
             /' psi_lm(r) = u_l(|r|) * Y_lm(r/|r|)')
 2010 format(/' augmentation functions for atom',i4,' (Rmt=',F6.3,')')
 2020 format(/' regular augmentation functions at r = Rmt' &
             /'  L    E(L)      u(r,E)     u''(r,E)   dE u(r,E)' &
             ,'  dE u''(r,E)  <u|dE u>  |dE u|^2')
 2030 format(I3,F8.3,1P,4E12.4,2E10.2)
 2070 format ('LO COEFFICIENT: l,A,B,C  ',i2,5X,6F12.5)
end subroutine auggen
end module     auggen_m


module     auglo_m; contains
subroutine auglo(stru, latom, rotloc, pos, il, Alm, Y, BK, coef, Nmat)
  use lolog,     only: iLO
  use loabc,     only: Alo
  use const,     only: DPk
  use param,     only: LOmax, Nrf
  use wplot,     only: Lmax7
  use structmod, only: struct_t

  !! procedure includes
  use Ylm_m

  implicit none

  type(struct_t), intent(in)    :: stru
  integer,        intent(in)    :: latom, Nmat
  integer,        intent(inout) :: il
  complex(DPk),   intent(inout) :: Alm((Lmax7+1)**2, Nrf)
  complex(DPk),   intent(in)    :: coef(Nmat)
  real(DPk),      intent(in)    :: &
       BK(3,Nmat), rotloc(3,3,stru%Nat), pos(3,stru%Nat)

  complex(DPk) :: PHS, PHSSUM(-LOMAX:LOMAX), Y((LOMAX+1)*(LOMAX+1))
  real(DPk)    :: RK(3), arg
  integer      :: l, lm, m, m1, irf, jlo, jatom, jneq

  jatom = stru%neq2at(latom)
  do L=0,LOmax
     do jlo=1,ilo(l,jatom)
        phssum = 0
        do jneq=1,stru%mult(jatom)
           do m1=-l,l
              il=il+1
              rk = matmul(rotloc(:,:,latom), BK(:,il))

              call YLM(rk, LOmax, Y)
              arg = dot_product(BK(:, il), pos(:, latom))

              phs=coef(il) * cmplx(cos(arg), sin(arg), DPk)

              LM = L**2
              do M=-L,L
                 LM = LM + 1
                 phssum(M) = phssum(M) + PHS * conjg( Y(LM) )
              ! write(0,'(i3, 2x, *(2e15.7, 2x))') lm, phs, y(lm)
              enddo
           enddo
        end do

        LM = L**2
        do M=-L,L
           LM=LM+1
           do irf=1,Nrf
              ALM(LM,irf)=ALM(LM,irf)+PHSSUM(M)*ALO(L,jlo,irf,JATOM)
           enddo
        enddo
     end do
  end do
  ! write(0,'(*(2e15.7, 2x))') sum(phssum), sum(alo)
end subroutine auglo
end module     auglo_m


module     augpw_m; contains
subroutine augpw(stru, latom, rotloc, pos, Npw, Alm, Y, BK, coef, Nmat)
  use atspdt,    only: P, DP
  use bessfu,    only: irad, fj, dfj
  use lolog,     only: lapw
  use param,     only: Nrf
  use wplot,     only: kconjg, Lmax7
  use const,     only: DPk
  use structmod, only: struct_t

  !! procedure includes
  use Ylm_m

  implicit none

  type(struct_t), intent(in)  :: stru
  integer,        intent(in)  :: Latom, Npw, Nmat
  real(DPk),      intent(in)  :: &
       BK(3,Nmat), rotloc(3,3,stru%Nat), pos(3,stru%Nat)
  complex(DPk),   intent(in)  :: COEF(Nmat)
  complex(DPk),   intent(out) :: Alm((LMAX7+1)**2, Nrf)
  complex(DPk),   intent(out) :: Y  ((LMAX7+1)**2)

  complex(DPk) :: PHS, PHSLM
  real(DPk)    :: RK(3), arg, al, bl
  integer      :: imt, jatom, iPW, lm, l, m

! The PW part of the augmentation coefficients A_lm,a and B_lm,a
! of a given eigen state at a given atom a in the unit cell
! ---------------------------------------------------------------------------
! Input:
! LATOM    -- the atom a
! NPW      -- current number of PW basis functions
! ROTLOC   -- the local rotation matrices T_a of each atom
!             [ stored as (T_a^-1)_ij = ROTLOC(i,j,a) ]
!
! Output:
! ALM  -- the PW part of the augmentation coefficients for each (l,m)
!
! Working arrays:
! Y        -- to hold Y_lm(...)
! ---------------------------------------------------------------------------
! X_l,m,a = Sum(K) c_K/sqrt(V) exp(i(K+k)R_a) Y(*)_lm(T_a^-1(K+k)) X_l,a(|K+k|)
!
! Here (*) stands for an optional complex conjugation on Y_lm(...)
! WIEN95 : (*) =
! WIEN97 : (*) = *
!
! R_a   : center of atom a
! T_a   : local rotation matrix of atom a
! X_l,a : PW augmentation coefficients for atom a
!
! Note: Letting  R_a = Q_a(R_0) + t_a  yields
!
!       exp(i(K+k)R_a) = exp(i(K+k)t_a + i[Q_a^-1(K+k)]R_0)
!
!       which is the expression used in LAPW1 and LAPW2
! ---------------------------------------------------------------------------

  jatom = stru%neq2at(latom)
  iMT   = iRad (jatom)

  !     << initialize ALM and BLM >>
  ALM = 0

  do IPW=1,NPW

     !       << Y*(*)_lm(T_a^-1(K+k)) >>
     RK = matmul(rotloc(:,:,latom), BK(:,ipw))

     call YLM(RK,LMAX7,Y)

     !:17[
     if(.not.KCONJG)then
        !         << WIEN95 convention : (*) =   >>
        do LM=1,(LMAX7+1)*(LMAX7+1)
           Y(LM) = conjg(Y(LM))
        end do
     end if
     !:17]

     !       << c_K/sqrt(V) * exp(i(K+k)R_a) >>
     arg = dot_product(BK(:, ipw), pos(:, latom))
     phs = coef(ipw) * cmplx(cos(arg), sin(arg), DPk)

     !       << A_l,a and B_l,a (without Rmt^2 factor) >>
     ! -----------------------------------------------------------------
     ! A_a,l = - [ d/dE d/dr u_l,a(Rmt,E)      j_l(Rmt*|K+k|) -
     !             d/dE      u_l,a(Rmt,E) d/dr j_l(Rmt*|K+k|) ]
     !
     ! B_a,l = - [           u_l,a(Rmt,E) d/dr j_l(Rmt*|K+k|) -
     !                  d/dr u_l,a(Rmt,E)      j_l(Rmt*|K+k|) ]
     ! -----------------------------------------------------------------
     LM = 0
     do L=0,LMAX7
        if (lapw(l,jatom)) then
           AL = DFJ(L,IPW,IMT)* P(L,2,JATOM) - &
                FJ(L,IPW,IMT)*DP(L,2,JATOM)
           BL =  FJ(L,IPW,IMT)*DP (L,1,JATOM) -  &
                DFJ(L,IPW,IMT)* P (L,1,JATOM)
        else
           AL = FJ(l,ipw,imt) / P(l,1,jatom) / stru%RMT(jatom)**2
           BL= 0.d0
        endif
        do M=-L,L
           LM = LM + 1
           PHSLM = PHS * conjg( Y(LM) )
           ALM(LM,1) = ALM(LM,1) + PHSLM * AL
           ALM(LM,2) = ALM(LM,2) + PHSLM * BL
        end do
     end do
  end do
end subroutine augpw
end module     augpw_m


module     latgen_m; contains
subroutine latgen(stru)
  use const,     only: DPk, TAU, SQ3
  use latt,      only: br1, br2, br3, br4
  use param,     only: unit_out
  use structmod, only: struct_t
  use clio,      only: croak

  !! procedure includes
  use gbass_m

  implicit none

  type(struct_t), intent(in) :: stru

  real(DPk)    :: alpha, beta, gamma, cosPhi, phi
  integer      :: i, j

! << Output (to module LATT) >>
! BR1(i,:) -- the real space lattice vectors a_i of the conventional unit cell
! BR2(i,:) -- the real space lattice vectors a_i of the primitive unit cell
! BR3(i,:) -- the reciprocal lattice vectors b_i of the conventional unit cell
! BR4(i,:) -- the reciprocal lattice vectors b_i of the primitive unit cell
!
! Here, the reciprocal lattice vectors are evaluated without a factor of 2pi !
!
! Caution: The lattice vectors setup here must precisely co-incide with
!          those used within LAPW2 !
!
  ALPHA = stru%alpha(1) * TAU/360
  BETA  = stru%alpha(2) * TAU/360
  GAMMA = stru%alpha(3) * TAU/360

  BR1=0; BR2=0

  if(STRU%LATTIC(1:1) == 'P')then
! -------------------------------------------------------------------------
!       << primitive lattice : P a b c alp bet gam >>
!
!       a_1 = a * (sin(bet)*sin(phi),sin(bet)*cos(phi),cos(bet))
!       a_2 = b * (        0        ,      sin(alp)   ,cos(alp))
!       a_3 = c * (        0        ,        0        ,  1     )
!       with
!       cos(phi) := ( cos(gam) - cos(alp)*cos(bet) ) / sin(alp)*sin(bet)
!
!       triclinic, monoclinic, orthorhombic, tetragonal, cubic
!
!       b_1 = ( 1/sin(bet)*sin(phi)                    ,     0     ,0) / a
!       b_2 = (-1/sin(alp)*tan(phi)                    , 1/sin(alp),0) / b
!       b_3 = ( 1/tan(alp)*tan(phi)-1/tan(bet)*sin(phi),-1/tan(alp),1) / c
! -------------------------------------------------------------------------
     COSPHI=(cos(GAMMA)-cos(ALPHA)*cos(BETA))/sin(ALPHA)/sin(BETA)
     PHI=acos(COSPHI)
!       << primitive unit cell >>
     BR2(1,1)=STRU%A(1)*sin(BETA)*sin(PHI)
     BR2(1,2)=STRU%A(1)*sin(BETA)*cos(PHI)
     BR2(1,3)=STRU%A(1)*cos(BETA)
     BR2(2,2)=STRU%A(2)*sin(ALPHA)
     BR2(2,3)=STRU%A(2)*cos(ALPHA)
     BR2(3,3)=STRU%A(3)
!       << conventional unit cell >>
     BR1 = BR2
  else if(STRU%LATTIC(1:1) == 'H') then
! -------------------------------------------------------------------------
!       << hexagonal lattice : H a * c * * * >>
!
!       a_1 = a * (sqrt(3)/2,-1/2,0)
!       a_2 = a * (     0   ,  1 ,0)
!       a_3 = c * (     0   ,  0 ,1)
!
!       this setting corresponds to the P a a c 90 90 120 setting
!
!       b_1 = ( 1/sqrt(3),  -1 ,   1 ) / a
!       b_2 = ( 1/sqrt(3),   1 ,   1 ) / a
!       b_3 = (-2/sqrt(3),   0 ,   1 ) / c
! -------------------------------------------------------------------------
!       << primitive unit cell >>
     BR2(1,1)= STRU%A(1)*SQ3/2
     BR2(1,2)=-STRU%A(1)/2
     BR2(2,2)= STRU%A(1)
     BR2(3,3)= STRU%A(3)
!       << conventional unit cell >>
     BR1 = BR2
  else if(STRU%LATTIC(1:1) == 'T'.or.STRU%LATTIC(1:1) == 'R') then
! -------------------------------------------------------------------------
!       << rhombohedral or trigonal lattice : R a * c * * * >>
!
!       a_1 = ( a/(2*sqrt(3)),-a/2,c/3)
!       a_2 = ( a/(2*sqrt(3)), a/2,c/3)
!       a_3 = (-a/   sqrt(3) ,  0 ,c/3)
!
!       Note: Although the trigonal lattice is treated as a primitive
!             lattice the lattice parameter correspond to the surrounding
!             (non-primitive) hexagonal unit cell.
!
!       b_1 = ( 1/(a*sqrt(3)), -1/b , 1/c )
!       b_2 = ( 1/(a*sqrt(3)),  1/b , 1/c )
!       b_3 = (-2/(a*sqrt(3)),   0  , 1/c )
! -------------------------------------------------------------------------
!       << primitive unit cell >>
     BR2(1,1)= STRU%A(1)/2/SQ3
     BR2(1,2)=-STRU%A(1)/2
     BR2(1,3)= STRU%A(3)/3
     BR2(2,1)= STRU%A(1)/2/SQ3
     BR2(2,2)= STRU%A(1)/2
     BR2(2,3)= STRU%A(3)/3
     BR2(3,1)=-STRU%A(1)/SQ3
     BR2(3,3)= STRU%A(3)/3
!       << conventional unit cell >>
     BR1 = BR2
  else if(STRU%LATTIC(1:1) == 'F') then
! -------------------------------------------------------------------------
!       << face-centered lattice : F a b c * * *
!
!       a_1 = ( 0 ,b/2,c/2)
!       a_2 = (a/2, 0 ,c/2)
!       a_3 = (a/2,b/2, 0 )
!
!       orthorhombic, cubic
!
!       b_1 = ( -1/a ,  1/b ,  1/c )
!       b_2 = (  1/a , -1/b ,  1/c )
!       b_3 = (  1/a ,  1/b , -1/c )
! -------------------------------------------------------------------------
!       << primitive unit cell >>
     BR2(1,2)=STRU%A(2)/2
     BR2(1,3)=STRU%A(3)/2
     BR2(2,1)=STRU%A(1)/2
     BR2(2,3)=STRU%A(3)/2
     BR2(3,1)=STRU%A(1)/2
     BR2(3,2)=STRU%A(2)/2
!       << conventional unit cell >>
     BR1(1,1)=STRU%A(1)
     BR1(2,2)=STRU%A(2)
     BR1(3,3)=STRU%A(3)
  else if(STRU%LATTIC(1:1) == 'B') then
! -------------------------------------------------------------------------
!       << body-centered lattice : B a b c * * *
!
!       a_1 = (-a/2, b/2, c/2)
!       a_2 = ( a/2,-b/2, c/2)
!       a_3 = ( a/2, b/2,-c/2)
!
!       orthorhombic, tetragonal, cubic
!
!       b_1 = (  0  , 1/b , 1/c )
!       b_2 = ( 1/a ,  0  , 1/c )
!       b_3 = ( 1/a , 1/b ,  0  )
! -------------------------------------------------------------------------
!       << primitive unit cell >>
     BR2(1,1)=-STRU%A(1)/2
     BR2(1,2)= STRU%A(2)/2
     BR2(1,3)= STRU%A(3)/2
     BR2(2,1)= STRU%A(1)/2
     BR2(2,2)=-STRU%A(2)/2
     BR2(2,3)= STRU%A(3)/2
     BR2(3,1)= STRU%A(1)/2
     BR2(3,2)= STRU%A(2)/2
     BR2(3,3)=-STRU%A(3)/2
!       << conventional unit cell >>
     BR1(1,1)=STRU%A(1)
     BR1(2,2)=STRU%A(2)
     BR1(3,3)=STRU%A(3)
  else if(stru%lattic(1:3)=='CXY' .and. stru%alpha(3)==90) then
! -------------------------------------------------------------------------
!       << base-centered (in the xy-plane) : CXY a b c alp bet 90 >>
!
!       Note: either alp or bet must be 90 degree
!
!       a_1 = (a*sin(bet)/2,-b*sin(alp)/2,a*cos(bet)/2-b*cos(alp)/2)
!       a_2 = (a*sin(bet)/2, b*sin(alp)/2,a*cos(bet)/2+b*cos(alp)/2)
!       a_3 = (    0       ,    0       ,            c            )
!
!       monoclinic, orthorhombic
!
!       b_1 = ( 1/(a*sin(bet)),-1/(b*sin(alp)), 0 )
!       b_2 = ( 1/(a*sin(bet)), 1/(b*sin(alp)), 0 )
!       b_3 = (-1/(c*tan(bet)),-1/(c*tan(alp)),1/c)
! -------------------------------------------------------------------------
     if(stru%alpha(1)/=90 .and. stru%alpha(2)/=90) &
          call croak('LATTIC NOT DEFINED: '//stru%lattic)
!       << primitive unit cell >>
     BR2(1,1)= STRU%A(1)/2*sin(BETA)
     BR2(1,2)=-STRU%A(2)/2*sin(ALPHA)
     BR2(1,3)= STRU%A(1)/2*cos(BETA) &
          - STRU%A(2)/2*cos(BETA)
     BR2(2,1)= STRU%A(1)/2*sin(BETA)
     BR2(2,2)= STRU%A(2)/2*sin(ALPHA)
     BR2(2,3)= STRU%A(1)/2*cos(BETA) &
          + STRU%A(2)/2*cos(ALPHA)
     BR2(3,3)= STRU%A(3)
!       << conventional unit cell >>
     BR1(1,1)=STRU%A(1)*sin(BETA)
     BR1(1,3)=STRU%A(1)*cos(BETA)
     BR1(2,2)=STRU%A(2)*sin(ALPHA)
     BR1(2,3)=STRU%A(2)*cos(ALPHA)
     BR1(3,3)=STRU%A(3)
  else if(stru%lattic(1:3)=='CYZ' .and. stru%alpha(1)==90) then
! -------------------------------------------------------------------------
!       << base-centered (in the yz-plane) : CYZ a b c 90 bet gam >>
!
!       Note: either bet or gam must be 90 degree
!
!       a_1 = (a*sin(bet)*sin(gam),a*sin(bet)*cos(gam),a*cos(bet))
!       a_2 = (        0        ,           b/2       ,  -c/2    )
!       a_3 = (        0        ,           b/2       ,   c/2    )
!
!       monoclinic, orthorhombic
!
!       b_1 = ( 1/(a*sin(bet)*sin(gam))               , 0 ,  0 )
!       b_2 = (-1/(b*tan(gam))+1/(c*tan(bet)*sin(gam)),1/b,-1/c)
!       b_3 = (-1/(b*tan(gam))-1/(c*tan(bet)*sin(gam)),1/b, 1/c)
! -------------------------------------------------------------------------
     if(stru%alpha(2)/=90 .and. stru%alpha(3)/=90) &
          call croak('LATTIC NOT DEFINED: '//stru%lattic)
!       << primitive unit cell >>
     BR2(1,1)= STRU%A(1)*sin(BETA)*sin(GAMMA)
     BR2(1,2)= STRU%A(1)*sin(BETA)*cos(GAMMA)
     BR2(1,3)= STRU%A(1)*cos(BETA)
     BR2(2,2)= STRU%A(2)/2
     BR2(2,3)=-STRU%A(3)/2
     BR2(3,2)= STRU%A(2)/2
     BR2(3,3)= STRU%A(3)/2
!       << conventional unit cell >>
     BR1(1,1)=STRU%A(1)*sin(BETA)*sin(GAMMA)
     BR1(1,2)=STRU%A(1)*sin(BETA)*cos(GAMMA)
     BR1(1,3)=STRU%A(1)*cos(BETA)
     BR1(2,2)=STRU%A(2)
     BR1(3,3)=STRU%A(3)
  else if(stru%lattic(1:3)=='CXZ' .and. stru%alpha(2)==90) then
! -------------------------------------------------------------------------
!       << base-centered (in the xz-plane) : CXZ a b c alp 90 gam >>
!
!       Note: either alp or gam must be 90 degree
!
!       a_1 = (a*sin(gam)/2,a*cos(gam)/2,  -c/2    )
!       a_2 = (    0       ,b*sin(alp)  ,b*cos(alp))
!       a_3 = (a*sin(gam)/2,a*cos(gam)/2,   c/2    )
!
!       monoclinic, orthorhombic
!
!       b_1 = ( 1/(a*sin(gam))         , 1/(c*tan(alp)),-1/c)
!       b_2 = (-1/(b*tan(gam)*sin(alp)), 1/(b*sin(alp)),  0 )
!       b_3 = ( 1/(a*sin(gam))         ,-1/(c*tan(alp)), 1/c)
! -------------------------------------------------------------------------
     if(stru%alpha(1)/=90 .and. stru%alpha(3)/=90) &
          call croak('LATTIC NOT DEFINED: '//stru%lattic)
!       << primitive unit cell >>
     BR2(1,1)= STRU%A(1)/2*sin(GAMMA)
     BR2(1,2)= STRU%A(1)/2*cos(GAMMA)
     BR2(1,3)=-STRU%A(3)/2
     BR2(2,2)= STRU%A(2)*sin(ALPHA)
     BR2(2,3)= STRU%A(2)*cos(ALPHA)
     BR2(3,1)= STRU%A(1)/2*sin(GAMMA)
     BR2(3,2)= STRU%A(1)/2*cos(GAMMA)
     BR2(3,3)= STRU%A(3)/2
!       << conventional unit cell >>
     BR1(1,1)=STRU%A(1)*sin(GAMMA)
     BR1(1,2)=STRU%A(1)*cos(GAMMA)
     BR1(2,2)=STRU%A(2)*sin(ALPHA)
     BR1(2,3)=STRU%A(2)*cos(ALPHA)
     BR1(3,3)=STRU%A(3)
  else
     stop 'LATTIC NOT DEFINED'
  end if

!     << find reciprocal lattice vectors (without a factor of 2pi) >>
  call GBASS(BR1,BR3,.false.)
  call GBASS(BR2,BR4,.false.)

  write(unit_out,1000) (I,(BR1(I,J),J=1,3),I=1,3)
  write(unit_out,1010) (I,(BR2(I,J),J=1,3),I=1,3)
  write(unit_out,1020) (I,(BR3(I,J),J=1,3),I=1,3)
  write(unit_out,1030) (I,(BR4(I,J),J=1,3),I=1,3)
1000 format(/' REAL SPACE LATTICE VECTORS a1, a2, a3 (in Bohr)' &
       /' -----------------------------------------------' &
       /' CONVENTIONAL UNIT CELL :'/(' a',I1,'     = ',3F12.7))
1010 format(/' PRIMITIVE UNIT CELL :'/(' a',I1,'     = ',3F12.7))
1020 format(/' RECIPROCAL LATTIC VECTORS b1, b2, b3 (in 1/Bohr)' &
       /' ------------------------------------------------' &
       /' CONVENTIONAL UNIT CELL :'/(' b',I1,'/2pi = ',3F12.7))
1030 format(/' PRIMITIVE UNIT CELL :'/(' b',I1,'/2pi = ',3F12.7))
end subroutine latgen
end module     latgen_m


module     wavsph_m; contains
subroutine wavsph(R, Bfac, iAt, iR, Psi, Y, stru)
  use work,      only: aug
  use const,     only: DPk
  use wplot,     only: kconjg, Lmax7
  use structmod, only: struct_t, radgrid

  !! procedure includes
  use Ylm_m

  implicit none

  real(DPk),      intent(in)  :: R(3)
  complex(DPk),   intent(in)  :: BFAC
  integer,        intent(in)  :: iAt, iR
  complex(DPk),   intent(out) :: Psi, Y((LMAX7+1)*(LMAX7+1))
  type(struct_t), intent(in)  :: stru

  complex(DPk) :: PHS, PHSLM
  integer      :: lm, jatom
  real(DPk)    :: RR, R1, R2, W1, W2

! evaluation of the wave function in the a muffin tin sphere around R+R_a
!
! psi(r) = e^ikR Sum(lm) w_lm,a(|r-R-R_a|) Y*(*)_lm(T_a^-1(r-R-R_a))
!
! Here (*) stands for an optional additional complex conjugation on Y*_lm(⋯)
! WIEN95 : (*) =     and hence   *(*) = *
! WIEN97 : (*) = *   and hence   *(*) =
! -----------------------------------------------------------------------
! Input:
! R    -- grid point in the local coordinate system of atom a
!         i.e the value s_a := T_a^-1(r-R-R_a)
! BFAC -- the Bloch factor e^ikR
! IAT  -- the atom a within the unit cell
! IR   -- the radial interval [r_i,r_i+1] the value |r-R-R_a| falls in
!
! MODULE WORK
! AUG(:,lm,a) -- the augmentation functions w_lm,a(r) on the radial mesh
!
! MODULE STRUCT  structural information
! MODULE RADGRD  radial grid information
!
! Output:
! PSI  -- the wave function psi(r)
!
! Working arrays:
! Y    -- to hold Y_lm(...)
! -----------------------------------------------------------------------

!     << Y*(*)_lm(T_a^-1(r-R-R_a)) for all lm >>
  call YLM(R,LMAX7,Y)
!:17[
  if(.not.KCONJG)then
!       << WIEN95 convention : (*) =   >>
     do lm=1,(LMAX7+1)**2
        Y(lm) = conjg(Y(lm))
     end do
  end if
!:17]
!
!     << prepare radial interpolation >>
! --------------------------------------------------------------------------
! r in [r1,r2] : w(r) = [ w(r1) ( r - r2) + w(r2) ( r1 - r ) ] / ( r2 - r1 )
! --------------------------------------------------------------------------
  jatom = stru%neq2at(iat)
  RR = sqrt( dot_product(R, R) )
  R1 = radgrid(stru, jatom, iR)
  R2 = radgrid(stru, jatom, iR+1)
  W1 = (R2-RR)/(R2-R1)
  W2 = (RR-R1)/(R2-R1)

!     << Sum(lm) ... >>
  PHS = 0
  do LM=1,(LMAX7+1)**2
     PHSLM = W1*AUG(IR,LM,IAT) + W2*AUG(IR+1,LM,IAT)
     PHS = PHS + PHSLM * Y(LM)
  end do

!     << psi(r) >>
  PSI = BFAC * PHS
end subroutine wavsph
end module     wavsph_m


module     bessel_m; contains
subroutine bessel(Npw, BK, NMT, rad, F, DF)
  use const, only: DPk
  use wplot, only: Lmax7

  !! procedure includes
  use sphbes_m
  use dvbes1_m

  implicit none

  ! Effective dimensions:
  !
  !   BK(3, Nmatk); rad(NMT); F(Lmax7+1, Nmatk, NMT),
  !   DF(Lmax7+1, Nmatk, NMT)
  !
  ! The actual array sizes in place of Nmatk and NMT may be different.

  integer,   intent(in)  :: Npw, NMT
  real(DPk), intent(in)  :: BK(:,:), rad(:)
  real(DPk), intent(out) :: F(:,:,:), DF(:,:,:)

  integer   :: ldim, ir, ik
  real(DPk) :: AK, arg, RMT

  ldim = size(F, 1)
  do iR=1,NMT
     RMT = rad(iR)
     do ik=1,Npw
        AK  = sqrt( dot_product(BK(:,ik), BK(:,ik)) )
        arg = AK * RMT
        call sphbes(ldim-1, arg, F(:, ik, iR))
        call dvbes1(F(:,ik,iR), DF(:, ik, iR), arg, ldim)
        DF(:, ik, iR) = AK*DF(:, ik, iR)
     end do
  end do
end subroutine bessel
end module     bessel_m


module     rotdef_m; contains
subroutine rotdef(stru, iop, pos)
  use const,     only: DPk
  use structmod, only: struct_t
  use clio,      only: croak
  use param,     only: unit_out
  use wplot,     only: mvatom
  use latt,      only: br2
  use sym2,      only: imat, trans, iord

  implicit none

  type(struct_t), intent(in)    :: stru
  integer,        intent(out)   :: iop(stru%Nneq*48)
  real(DPk),      intent(inout) :: pos(3, stru%Nat)

  real(DPk), parameter :: half = 0.5_DPk, tol = 1e-4_DPk

  integer   :: index, index1, jatom, jOp, m
  real(DPk) :: pos0(3), R(3), D(3)

!!! << Input >>
!!! stru     -- struct
!!!
!!! from MODULE SYM2 -- symmetry operations
!!! IMAT,TRANS -- symmetry operations (in primitive fractional coordinates)
!!! IORD     -- number of symmetry operations
!!! -------------------------------------------------------------------------
!!! {Q|t} : y_i = Sum(j) Q_ij x_j + t_i  with  Q_ij = IMAT(j,i) and t_i = TRANS(i)
!!! -------------------------------------------------------------------------
!!!
!!! << Output >>
!!! POS(:,i) -- the proper symmetry generated position of the i-th atom
!!! IOP(i)   -- the symmetry operation which generates the i-th atom from
!!!             the first atom of the corresponding set of inequivalent atoms
!!!
!!! Caution: The symmetry operations selected here to generate the individual
!!!          atoms must precisely co-incide with those selected in LAPW2 !

!!!     << find the generating symmetry operations >>
  write(unit_out,2000)

  index=0
  nneq: do jatom=1,stru%Nneq
     index1 = index+1
!!!       << store r_1.atom before updating it >>
     pos0 = pos(:,index1)

     mult: do M=1,stru%mult(jatom)
        index = index+1
        symop: do jop=1,iord

!!!           << find {Q|t}(r_1.atom) >>
           R = matmul(transpose(imat(:,:,jOp)), pos0 + trans(:,jOp))

!!!           << check the difference modulo lattice translations >>
           D = abs(mod(abs(pos(:,index) - R) + HALF, 1._DPk) - HALF)

           if (D(1).lt.TOL.and.D(2).lt.TOL.and.D(3).lt.TOL) then
              IOP(INDEX) = JOP
              !:17[
              if(.not. mvatom)then
                 R = pos(:, index)
              endif
              !:17]
              pos(:, index) = R

!!!             << print updated position in Cartesian coordinates >>
              R = matmul(BR2, pos(:, index))

              write(unit_out,2010) jatom, IOP(index), R
              cycle mult
           end if
        end do symop

!!!         << something is wrong here >>
        write(unit_out,1000) index1, index
        write(unit_out,1010) index1, pos(:, index1)
        write(unit_out,1010) index , pos(:, index )
        call croak('error in ROTDEF')
     end do mult
  end do nneq

1000 format(///3X,'ERROR IN ROTDEF:'/ &
       'NO SYMMETRY OPERATION FOUND TO MAP ATOM',I3, &
       ' ONTO ATOM',I3)
1010 format('ATOM',I3,' : POS = ',3F13.7)

2000 format(/' SYMMETRY ADJUSTED POSITIONS OF THE BASIS ATOMS' &
       /' ----------------------------------------------' &
       /' atom  symm.     x [Bohr]     y [Bohr]     z [Bohr]')
2010 format(2(I5,1X),3F13.7)
end subroutine rotdef
end module     rotdef_m


module     findmt_m; contains
subroutine findmt(P, atms, stru, pos, iAt, iLat, iR, R)
  use const,     only: DPk
  use latt,      only: br2
  use structmod, only: struct_t

  implicit none

  type(struct_t), intent(in)  :: stru
  real(DPk),      intent(in)  :: P(3), atms(3, stru%Nneq), pos(3, stru%Nat)
  integer,        intent(out) :: iAt, iLat(3), iR
  real(DPk),      intent(out) :: R(3)

  integer   :: i, j, jatom, jx, jy, jz, ilow, iup
  real(DPk) :: RR, R2, RMT2, T

!
! checks whether a given real space point r falls into a muffin tin
! sphere around R + R_a with a running over all atoms in the unit cell
! --------------------------------------------------------------------
! Input:
! P(:)  -- the real space point in primitive fractional coordinates
! NAT   -- total number of atoms per unit cell
! ATMS  -- the size of the smallest primitive unit cells surrounding
!          the muffin tin spheres of each group of sym.-eq. atoms
!
! module LATT    data on the Bravais lattice
!
! Output:
! IAT   -- the atom a the muffin tin sphere belongs to
!          or zero if the given point is in the interstitial region
! ILAT  -- the lattice displacement R = Sum(i) N_i a_i in primitive
!          fractional coordinates [ undefined if IAT = 0 ]
! IR    -- radial mesh index i such that |r - R - R_a| in [r_i,r_i+1]
!          [ undefined if IAT = 0 ]
! R(:)  -- if in interstitial:
!          the absolute position r in Cartesian coordinates
!          if in muffin tin spheres:
!          the relative positition r - R - R_a in Cartesian corrdinates
!
!          Here a_i are the primitive real space lattice vectors !
! --------------------------------------------------------------------

  dimension T(3),ILOW(3),IUP(3)

  do iat=1,stru%Nat
     jatom = stru%neq2at(iat)
     RMT2 = stru%RMT(jatom)**2
!
!       << setup search area for atomic sphere at R0+R >>
! --------------------------------------------------------------------
! Letting R = Sum(i) N_i a_i the search area for R is given by
! all integers N_i in [ x_i - R0_i - s_i , x_i - R0_i + s_i ]
!
! to account for numerical noise let s_i = s_i + 0.0001 here!
! --------------------------------------------------------------------
     do I=1,3
        T   (I) = P(I) - pos(i,iat)
        ILOW(I) = nint( T(I) - ATMS(I,JATOM) + 0.4999_DPk )
        IUP (I) = nint( T(I) + ATMS(I,JATOM) - 0.4999_DPk )
     end do

!       << check all relevent atomic sphere displacements >>
     do JZ=ILOW(3),IUP(3)
        do JY=ILOW(2),IUP(2)
           do JX=ILOW(1),IUP(1)

!             << load x - R0 - R in Cartesian coordinates >>
              do J=1,3
                 R(J) = (T(1)-JX)*BR2(1,J) + (T(2)-JY)*BR2(2,J) &
                      + (T(3)-JZ)*BR2(3,J)
              end do

              R2 = dot_product(R,R)
              if( R2.lt.RMT2 ) then
!               << in muffin tin sphere IAT >>
                 RR = sqrt(R2)
                 ILAT(1) = JX
                 ILAT(2) = JY
                 ILAT(3) = JZ
                 IR = min(1 + int(log(max(RR/stru%r0(JATOM),1.0D0)) &
                      &           / stru%dx(JATOM)), &
                      &   stru%Npt(JATOM))
                 return
              end if
           end do
        end do
     end do
  end do

!     << in interstitial >>
  IAT = 0
  do J=1,3
     R(J) = P(1)*BR2(1,J) + P(2)*BR2(2,J) + P(3)*BR2(3,J)
  end do

  return
end subroutine findmt
end module     findmt_m


module     locdef_m; contains
subroutine locdef(rot0, imat, rot)
  use latt,  only: br2, br4
  use const, only: DPk
  !:17[
  use wplot, only: addloc, userot
  !:17]

  implicit none

  real(DPk), intent(in)  :: rot0(3,3)
  integer,   intent(in)  :: iMat(3,3)
  real(DPk), intent(out) :: rot (3,3)

  ! Input:
  ! ROT0 : local reference rotation matrix (Cartesian coordinates)
  !        x'_i = Sum(j) (T^-1)_ij x_j  with  T_ji = (T^-1)_ij = ROT0(i,j)
  ! IMAT : symmetry operation {Q|t} (primitive fractional coordinates)
  !        y_m = Sum(n) Q_mn x_n + t_m  with  Q_mn = IMAT(n,m)
  !
  ! from module LATT
  ! BR2  : primitive real space lattice vector am = BR2(m,:)
  ! BR4  : primitive reciprocal lattice vector bn = BR4(n,:) [without 2pi]
  !
  ! Output:
  ! ROT  : symmetry adapted local rotation matrix (Cartesian coordinates)
  !        x'_i = Sum(j) (R^-1)_ij x_j  with  R_ji = (R^-1)_ij = ROT(i,j)
  !        [ ROT and ROT0 may share memory ! ]
  !
  ! Algorithm:
  ! R := Q o T  with  Q_ik = Sum(m,n) am_i Q_mn bn_k

  real(DPk) :: Amat(3,3), Bmat(3,3)
  integer   :: k

  !:17[
  if(.not. userot)then
     ! don't apply any local rotations, i.e. R := E
     rot = 0

     do k=1,3
        rot(k,k) = 1
     end do
     return
  endif
  !:17]

  ! A_mk := Sum(n) Q_mn bn_k
  Amat = matmul(transpose(imat), BR4)

  !:17[
  if(.not. addloc)&
       ! ignore local rotation matrix T from input, i.e. R = Q
       ! R_ik = Q_ik = Sum(m) am_i A_mk
       rot = matmul(transpose(Amat), BR2)
  !:17]

  ! B_mj = Sum(k) A_mk T_kj
  Bmat = matmul(Amat, transpose(rot0))

  ! R_ij = Sum(m) am_i B_mj
  rot = matmul(transpose(Bmat), BR2)
end subroutine locdef
end module     locdef_m


module     trans_m; contains
subroutine trans(pos)
  use const,     only: DPk
  use latt,      only: br1, br2, br3, br4
  use wplot,     only: Nsym
  use sym2,      only: rtrans=>trans, imat

  implicit none

  real(DPk), intent(inout) :: pos(:, :)

  real(DPk) :: F(3), S(3,3), T(3,3), Q(3,3)
  integer   :: i, k, n

! transforms real space vectors x and symmetry operations {Q|t}
! from conventional into primitive fractional coordinates
! --------------------------------------------------------------------
! Input:
! NPOS     -- the number of real space vectors x to transform
! NSYM     -- the number of symmetry operations {Q|t} to transform
!
! module LATT
! BR1(i,:) -- conventional real space lattice vectors a_i
! BR2(i,:) -- primitive    real space lattice vectors a_i
! BR3(i,:) -- conventional reciprocal lattice vectors b_i (without 2pi)
! BR4(i,:) -- primitive    reciprocal lattice vectors b_i (without 2pi)
!
! Input/Output:
! POS(:,n)    -- the n-th real space vector x
!                x = Sum(i=1,3) x_i a_i  with  x_i = POS(i,n)
! IMAT(:,:,n) -- the n-th symmetry operation {Q,t} :
! RTRANS(:,n)       y_i = Sum(j) Q_ij x_i + t_i  with  Q_ij = IMAT(j,i,n)
!                                             and   t_i  = RTRANSK(  i,n)
!
! Algorithm:
! real space vectors from conventional to primitive:
! p_k = Sum(i) T(k,i) c_i  with  T(k,i) = Sum(j) BR4(k,j) BR1(i,j)
!
! real space vectors from primitive to conventional:
! c_i = Sum(k) S(i,k) p_k  with  S(i,k) = Sum(j) BR3(i,j) BR2(k,j)
!
! symmetry operations from conventional to primitive
! Q_kk' = Sum(ii') T(k,i) Q_ii' S(i',k')  and  t_k = Sum(i) T(k,i) t_i
! --------------------------------------------------------------------

!     << set up the transformation matrices >>
  T = matmul(BR4, transpose(BR1))
  S = matmul(BR3, transpose(BR2))
  write(0,*) 'trans()'
  write(0,'("T1=", 3F7.3)') T
  write(0,'("S1=", 3F7.3)') S

  do K=1,3
     do I=1,3
        T(K,I) = BR4(K,1)*BR1(I,1) + BR4(K,2)*BR1(I,2) &
             + BR4(K,3)*BR1(I,3)
        S(I,K) = BR3(I,1)*BR2(K,1) + BR3(I,2)*BR2(K,2) &
             + BR3(I,3)*BR2(K,3)
     end do
  end do

  write(0,'("T2=", 3F7.3)') T
  write(0,'("S2=", 3F7.3)') S

!     << transform the real space vectors >>
  do n=1,size(pos,2)
     pos(:, N) = matmul(T, pos(:, n))
  end do

!     << transform the symmetry operations >>
  do N=1,NSYM
     F = matmul(T, rtrans(:, N))
     Q = matmul(T, transpose(imat(:,:,N)))
     write(0,'("F1=", 3F7.3)') F
     write(0,'("Q1=", 3F7.3)') Q
     do K=1,3
        F(K) = T(K,1)*RTRANS(1,N) + T(K,2)*RTRANS(2,N) &
             + T(K,3)*RTRANS(3,N)

        do I=1,3
           Q(K,I) = T(K,1)*IMAT(I,1,N) + T(K,2)*IMAT(I,2,N) &
                & + T(K,3)*IMAT(I,3,N)
        end do
     end do
     write(0,'("F2=", 3F7.3)') F
     write(0,'("Q2=", 3F7.3)') Q

     rtrans(:, N) = F
     imat (:,:,N) = nint(matmul(Q, S))
  end do
end subroutine trans
end module     trans_m


!!/---
!! Local Variables:
!! mode: f90
!! End:
!!\---
!!
!! Time-stamp: <2016-11-29 16:55:33 assman@faepop71.tu-graz.ac.at>
