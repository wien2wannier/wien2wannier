!!! wien2wannier/SRC_w2w/util_w2w.F
!!!
!!!    Collection of routines for wien2wannier
!!!
!!! Copyright 2010-2012 Philipp Wissgott
!!!           2013-2015 Elias Assmann

!!/=== Things that are exported here =============
!!
!! const: Constants mathematical and configurational
!!
!!    BUFSZ, DPk, IPk, R4, R8, C16
!!    PI, ORTHO_TEST, SQ3, Ryd_eV, Ang_Bohr
!!
!!
!! reallocate: Reallocation of arrays of various types
!!
!!    realloc()
!!
!!
!! util: Assorted miscellanea
!!
!!    string(), uppercase(), lowercase(), inverse3x3()
!!    ptime(), type(ptimer), ptick(), ptock(), ptot(),
!!    line_count(), paropen(), newunit()
!!
!!
!! clio: Helper procedures for command line interface
!!
!!    croak(), carp(), cluck(), fetcharg(), type(argstr)
!!
!!
!! structmod: ‘struct_t’ type and associated procedures
!!
!!    type(struct_t), struct_read()
!!
!!
!! inwfmod: ‘inwf_t’ type and associated procedures
!!
!!    type(inwf_t), inwf_read(), type(proj)
!!
!!
!! kpoints: Reading ‘klist’ files
!!
!!    get_kmesh_klist(), get_kmesh_band(),
!!    count_kmesh_klist(), count_kmesh_band()
!!
!!
!! wannier90: Reading ‘chk’ files
!!
!!    type(chk_t), chk_read()
!!
!!
!! wien2k: Procedures nicked from Wien2k
!!
!!    errflg(), errclr(), gtfnam()
!!
!!\===============================================


!---------------  Constants mathematical and configurational  ---------------
module const
  use iso_fortran_env, only: error_unit

  implicit none
  public

  integer, parameter   :: BUFSZ = 256
  !! Default kinds
  integer, parameter   :: DPk = kind(1.0d0) ! inherited from W90 2.0.0
  integer, parameter   :: IPk = selected_int_kind(8)
  !! Kinds for WIEN-compatibility
  real*4,     private  :: four_real
  real*8,     private  :: eight_real
  complex*16, private  :: sixteen_complex
  integer, parameter   :: R4 = kind(four_real)
  integer, parameter   :: R8 = kind(eight_real)
  integer, parameter   :: C16= kind(sixteen_complex)

  real(DPk), parameter :: PI         = acos(-1._DPk)
  real(DPk), parameter :: ORTHO_TEST = 1.e-12_DPk
  real(DPk), parameter :: SQ3        = sqrt(3._DPk)
  real(DPk), parameter :: Ryd_eV     = 13.605698_DPk ! value for WIEN compat.
  real(DPk), parameter :: Ang_Bohr   = 1.8897261246_DPk
end module const


!----------------  Reallocation of various arrays           -----------------
module reallocate
  use const

  implicit none

  private
  public :: realloc

  interface realloc
     ! Fortran makes this pretty cumbersome.  We need to provide one
     ! specific subroutine for each type and rank
     module procedure realloc_r_d1, realloc_r_d2, realloc_r_d3
     module procedure realloc_c_d1, realloc_c_d2, realloc_c_d3
     module procedure realloc_i_d1, realloc_i_d2, realloc_i_d3
#ifdef HAVE_POINTER_ALLOCATABLE_GENERICS
     ! In Fortran 2008, generic procedures can distinguish between
     ! ‘pointer’ and ‘allocatable’ variables, but older compilers do
     ! not support this yet, so we only support ‘allocatable’ there.
     module procedure realloc_r_p1, realloc_r_p2, realloc_r_p3
     module procedure realloc_c_p1, realloc_c_p2, realloc_c_p3
     module procedure realloc_i_p1, realloc_i_p2, realloc_i_p3
#endif
  end interface realloc
contains
#ifdef HAVE_POINTER_ALLOCATABLE_GENERICS
!!! Procedures for ‘pointer’ variables, which use pointer assignment
  subroutine realloc_r_p1(tf, nd)
    integer,           intent(in)  :: nd
    ! pointer INTENT refers to association staus
    real(r8), pointer, intent(out) :: tf(:)

    real(r8), pointer :: hilfsfeld(:)
    integer           :: min1

    allocate(hilfsfeld(nd))
    !     nur 1 mal kopieren reicht
    !     auch fuer mehrdimensionale Felder schaut die Zuweisung gleich aus
    min1=min(nd,size(tf,1))
    hilfsfeld(1:min1)=tf(1:min1)
    deallocate(tf)
    !     der Zeiger wird nur auf das neue Feld umgebogen, nicht neu alloziert
    tf=>hilfsfeld
  end subroutine realloc_r_p1

  subroutine realloc_r_p2(tf, nd1, nd2)
    integer,           intent(in)  :: nd1, nd2
    real(r8), pointer, intent(out) :: tf(:,:)

    real(r8), pointer :: hilfsfeld(:,:)
    integer           :: min1, min2

    allocate(hilfsfeld(nd1,nd2))

    min1=min(nd1,size(tf,1))
    min2=min(nd2,size(tf,2))
    hilfsfeld(1:min1,1:min2)=tf(1:min1,1:min2)

    deallocate(tf)
    tf=>hilfsfeld
  end subroutine realloc_r_p2

  subroutine realloc_r_p3(tf, nd1, nd2, nd3)
    integer,           intent(in)  :: nd1, nd2, nd3
    real(r8), pointer, intent(out) :: tf(:,:,:)

    real(r8), pointer :: hilfsfeld(:,:,:)
    integer           :: min1, min2, min3

    allocate(hilfsfeld(nd1,nd2,nd3))

    min1=min(nd1,size(tf,1))
    min2=min(nd2,size(tf,2))
    min3=min(nd3,size(tf,3))
    hilfsfeld(1:min1,1:min2,1:min3)=tf(1:min1,1:min2,1:min3)

    deallocate(tf)
    tf=>hilfsfeld
  end subroutine realloc_r_p3

  subroutine realloc_c_p1(tf, nd)
    integer,               intent(in)  :: nd
    complex(c16), pointer, intent(out) :: tf(:)

    complex(c16), pointer :: hilfsfeld(:)
    integer               :: min1

    allocate(hilfsfeld(nd))

    min1=min(nd,size(tf,1))
    hilfsfeld(1:min1)=tf(1:min1)

    deallocate(tf)
    tf=>hilfsfeld
  end subroutine realloc_c_p1

  subroutine realloc_c_p2(tf, nd1, nd2)
    integer,               intent(in)  :: nd1, nd2
    complex(c16), pointer, intent(out) :: tf(:,:)

    complex(c16), pointer :: hilfsfeld(:,:)
    integer             :: min1, min2

    allocate(hilfsfeld(nd1,nd2))

    min1=min(nd1,size(tf,1))
    min2=min(nd2,size(tf,2))
    hilfsfeld(1:min1,1:min2)=tf(1:min1,1:min2)

    deallocate(tf)
    tf=>hilfsfeld
  end subroutine realloc_c_p2

  subroutine realloc_c_p3(tf, nd1, nd2, nd3)
    integer,               intent(in)  :: nd1, nd2, nd3
    complex(c16), pointer, intent(out) :: tf(:,:,:)

    complex(c16), pointer :: hilfsfeld(:,:,:)
    integer               :: min1, min2, min3

    allocate(hilfsfeld(nd1,nd2,nd3))

    min1=min(nd1,size(tf,1))
    min2=min(nd2,size(tf,2))
    min3=min(nd3,size(tf,3))
    hilfsfeld(1:min1,1:min2,1:min3)=tf(1:min1,1:min2,1:min3)

    deallocate(tf)
    tf=>hilfsfeld
  end subroutine realloc_c_p3

  subroutine realloc_i_p1(tf, nd)
    integer,          intent(in)  :: nd
    integer, pointer, intent(out) :: tf(:)

    integer, pointer :: hilfsfeld(:)
    integer          :: min1

    allocate(hilfsfeld(nd))

    min1=min(nd,size(tf,1))
    hilfsfeld(1:min1)=tf(1:min1)

    deallocate(tf)
    tf=>hilfsfeld
  end subroutine realloc_i_p1

  subroutine realloc_i_p2(tf, nd1, nd2)
    integer,          intent(in)  :: nd1, nd2
    integer, pointer, intent(out) :: tf(:,:)

    integer, pointer :: hilfsfeld(:,:)
    integer          :: min1, min2

    allocate(hilfsfeld(nd1,nd2))

    min1=min(nd1,size(tf,1))
    min2=min(nd2,size(tf,2))
    hilfsfeld(1:min1,1:min2)=tf(1:min1,1:min2)

    deallocate(tf)
    tf=>hilfsfeld
  end subroutine realloc_i_p2

  subroutine realloc_i_p3(tf, nd1, nd2, nd3)
    integer,          intent(in)  :: nd1, nd2, nd3
    integer, pointer, intent(out) :: tf(:,:,:)

    integer, pointer :: hilfsfeld(:,:,:)
    integer          :: min1, min2, min3

    allocate(hilfsfeld(nd1,nd2,nd3))

    min1=min(nd1,size(tf,1))
    min2=min(nd2,size(tf,2))
    min3=min(nd3,size(tf,3))
    hilfsfeld(1:min1,1:min2,1:min3)=tf(1:min1,1:min2,1:min3)

    deallocate(tf)
    tf=>hilfsfeld
  end subroutine realloc_i_p3
#endif

!!! Procedures for ‘allocatable’ variables, which use move_alloc()
  subroutine realloc_r_d1(tf, nd)
    integer,               intent(in)    :: nd
    real(r8), allocatable, intent(inout) :: tf(:)
                         
    real(r8), allocatable :: hilfsfeld(:)
    integer               :: min1

    allocate(hilfsfeld(nd))
    min1=min(nd,size(tf,1))
    hilfsfeld(1:min1)=tf(1:min1)

    call move_alloc(hilfsfeld, tf) ! hilfsfeld is automatically deallocated
  end subroutine realloc_r_d1

  subroutine realloc_r_d2(tf, nd1, nd2)
    integer,               intent(in)    :: nd1, nd2
    real(r8), allocatable, intent(inout) :: tf(:,:)
                         
    real(r8), allocatable :: hilfsfeld(:,:)
    integer               :: min1, min2

    allocate(hilfsfeld(nd1,nd2))

    min1=min(nd1,size(tf,1))
    min2=min(nd2,size(tf,2))
    hilfsfeld(1:min1,1:min2)=tf(1:min1,1:min2)

    call move_alloc(hilfsfeld, tf)
  end subroutine realloc_r_d2

  subroutine realloc_r_d3(tf, nd1, nd2, nd3)
    integer,               intent(in)    :: nd1, nd2, nd3
    real(r8), allocatable, intent(inout) :: tf(:,:,:)
                         
    real(r8), allocatable :: hilfsfeld(:,:,:)
    integer               :: min1, min2, min3

    allocate(hilfsfeld(nd1,nd2,nd3))

    min1=min(nd1,size(tf,1))
    min2=min(nd2,size(tf,2))
    min3=min(nd3,size(tf,3))
    hilfsfeld(1:min1,1:min2,1:min3)=tf(1:min1,1:min2,1:min3)

    call move_alloc(hilfsfeld, tf)
  end subroutine realloc_r_d3

  subroutine realloc_c_d1(tf, nd)
    integer,                   intent(in)    :: nd
    complex(c16), allocatable, intent(inout) :: tf(:)
                             
    complex(c16), allocatable :: hilfsfeld(:)
    integer                   :: min1

    allocate(hilfsfeld(nd))

    min1=min(nd,size(tf,1))
    hilfsfeld(1:min1)=tf(1:min1)

    call move_alloc(hilfsfeld, tf)
  end subroutine realloc_c_d1

  subroutine realloc_c_d2(tf, nd1, nd2)
    integer,                   intent(in)    :: nd1, nd2
    complex(c16), allocatable, intent(inout) :: tf(:,:)
                             
    complex(c16), allocatable :: hilfsfeld(:,:)
    integer                  :: min1, min2

    allocate(hilfsfeld(nd1,nd2))

    min1=min(nd1,size(tf,1))
    min2=min(nd2,size(tf,2))
    hilfsfeld(1:min1,1:min2)=tf(1:min1,1:min2)

    call move_alloc(hilfsfeld, tf)
  end subroutine realloc_c_d2

  subroutine realloc_c_d3(tf, nd1, nd2, nd3)
    integer,                   intent(in)    :: nd1, nd2, nd3
    complex(c16), allocatable, intent(inout) :: tf(:,:,:)
                             
    complex(c16), allocatable :: hilfsfeld(:,:,:)
    integer                   :: min1, min2, min3

    allocate(hilfsfeld(nd1,nd2,nd3))

    min1=min(nd1,size(tf,1))
    min2=min(nd2,size(tf,2))
    min3=min(nd3,size(tf,3))
    hilfsfeld(1:min1,1:min2,1:min3)=tf(1:min1,1:min2,1:min3)

    call move_alloc(hilfsfeld, tf)
  end subroutine realloc_c_d3

  subroutine realloc_i_d1(tf, nd)
    integer,              intent(in)    :: nd
    integer, allocatable, intent(inout) :: tf(:)
                        
    integer, allocatable :: hilfsfeld(:)
    integer              :: min1

    allocate(hilfsfeld(nd))

    min1=min(nd,size(tf,1))
    hilfsfeld(1:min1)=tf(1:min1)

    call move_alloc(hilfsfeld, tf)
  end subroutine realloc_i_d1

  subroutine realloc_i_d2(tf, nd1, nd2)
    integer,              intent(in)    :: nd1, nd2
    integer, allocatable, intent(inout) :: tf(:,:)
                        
    integer, allocatable  :: hilfsfeld(:,:)
    integer               :: min1, min2

    allocate(hilfsfeld(nd1,nd2))

    min1=min(nd1,size(tf,1))
    min2=min(nd2,size(tf,2))
    hilfsfeld(1:min1,1:min2)=tf(1:min1,1:min2)

    call move_alloc(hilfsfeld, tf)
  end subroutine realloc_i_d2

  subroutine realloc_i_d3(tf, nd1, nd2, nd3)
    integer,              intent(in)    :: nd1, nd2, nd3
    integer, allocatable, intent(inout) :: tf(:,:,:)

    integer, allocatable :: hilfsfeld(:,:,:)
    integer              :: min1, min2, min3

    allocate(hilfsfeld(nd1,nd2,nd3))

    min1=min(nd1,size(tf,1))
    min2=min(nd2,size(tf,2))
    min3=min(nd3,size(tf,3))
    hilfsfeld(1:min1,1:min2,1:min3)=tf(1:min1,1:min2,1:min3)

    call move_alloc(hilfsfeld, tf)
  end subroutine realloc_i_d3
end module reallocate


!---------------  Assorted miscellanea                        ---------------
MODULE util
  use const, only: DPk, BUFSZ, ERROR_unit
  implicit none

  private
  public :: string, inverse3x3, ptime, ptimer, ptick, ptock, ptot, &
       &    uppercase, lowercase, line_count, paropen, newunit

  interface string
     module procedure int2str, ints2str, real2str, reals2str
  end interface string

  interface ptime
!!! Ptime() takes and outputs CPU and wall time measurements.  Typical
!!! usage is
!!!
!!!	ptime([UNIT=…])
!!!
!!! to set output unit and initialize the timers; then
!!!
!!!	… do stuff …
!!!
!!!	call ptime("doing stuff")
!!!
!!!	… do other stuff …
!!!
!!!	call ptime("other stuff")
     module procedure ptime_setunit, ptime_print, ptick
  end interface

  type ptimer
     real(DPk)            :: cpu1 =0, cputot =0
     integer              :: wall1=0, walltot=0, countrate=0
     character(len=BUFSZ) :: &
          fmt = "('Times for ',A,T33,'(sec):',F8.3,' wall;',F9.3,' CPU')"
     integer              :: unit=ERROR_UNIT
     ! flush after every write? -- initialized to F here because we do
     ! not know if the target unit will support it
     logical              :: flush=.false.
  end type ptimer

  type(ptimer), target, save :: default_timer

contains
  pure character(len=15) function int2str(n)
    integer, intent(in) :: n
    write(int2str,"(I0)") n
  end function int2str

  pure character(len=BUFSZ) function ints2str(na)
    integer, intent(in) :: na(:)

    character(BUFSZ) :: fmt

    write(fmt, '( "( ''[ ''", I0, "(I0,1X)", "'' ]'' )" )') size(na)
    write(ints2str, fmt) na
  end function ints2str

  pure character(len=15) function real2str(r)
    real(DPk), intent(in) :: r
    write(real2str,"(E15.8)") r
  end function real2str

  pure character(len=BUFSZ) function reals2str(ra)
    real(DPk), intent(in) :: ra(:)

    character(BUFSZ) :: fmt

    write(fmt, '( "( ''[''", I0, "(F0.4,1X)", "'' ]'' )" )') size(ra)
    write(reals2str, fmt) ra
  end function reals2str

  pure subroutine inverse3x3(a, ainv)
    !inverse of th 3x3 matrix A
    implicit none

    real(DPk), intent(in)  :: a(3,3)
    real(DPk), intent(out) :: ainv(3,3)
    real(DPk) :: det

    det = a(1,1)*a(2,2)*a(3,3) + a(1,2)*a(2,3)*a(3,1) &
         +a(1,3)*a(2,1)*a(3,2) - a(3,1)*a(2,2)*a(1,3) &
         -a(1,1)*a(3,2)*a(2,3) - a(2,1)*a(1,2)*a(3,3)

    ainv(1,1) = (  a(2,2)*a(3,3) - a(2,3)*a(3,2) ) / det
    ainv(2,1) = (- a(2,1)*a(3,3) + a(2,3)*a(3,1) ) / det
    ainv(3,1) = (  a(2,1)*a(3,2) - a(2,2)*a(3,1) ) / det
    ainv(1,2) = (- a(1,2)*a(3,3) + a(1,3)*a(3,2) ) / det
    ainv(2,2) = (  a(1,1)*a(3,3) - a(1,3)*a(3,1) ) / det
    ainv(3,2) = (- a(1,1)*a(3,2) + a(1,2)*a(3,1) ) / det
    ainv(1,3) = (  a(1,2)*a(2,3) - a(1,3)*a(2,2) ) / det
    ainv(2,3) = (- a(1,1)*a(2,3) + a(1,3)*a(2,1) ) / det
    ainv(3,3) = (  a(1,1)*a(2,2) - a(1,2)*a(2,1) ) / det
  end subroutine inverse3x3

!!! ptime stuff
  subroutine ptime_setunit(unit, timer, flush)
    integer,              intent(in)              :: unit
    type(ptimer), target, intent(inout), optional :: timer
    logical,              intent(in),    optional :: flush

    integer               :: ios
    type(ptimer), pointer :: t
    t                     => default_timer
    if (present(timer)) t => timer

    t%unit = unit
    if (present(flush)) then
       t%flush = flush
    else                ! by default, we flush if possible
       flush(UNIT=unit, IOSTAT=ios)

       t%flush = ios==0 ! ios>0 signals error, ios<0 means flush not supported
    end if
       
    call ptick(t)
  end subroutine ptime_setunit

  subroutine ptick(timer)
    type(ptimer), target, intent(inout), optional :: timer
    
    type(ptimer), pointer :: t
    t                      => default_timer
    if (present(timer)) t  => timer

    call cpu_time    (t% cpu1)
    call system_clock(t%wall1, t%countrate)
  end subroutine ptick

  subroutine ptock(timer)
    type(ptimer), target, intent(inout), optional :: timer
    
    real(DPk)             :: cpu2
    integer               :: wall2
    type(ptimer), pointer :: t
    t                      => default_timer
    if (present(timer)) t  => timer

    call cpu_time    ( cpu2)
    call system_clock(wall2)

    t% cputot = t% cputot +  cpu2 - t% cpu1
    t%walltot = t%walltot + wall2 - t%wall1

    t% cpu1 =  cpu2
    t%wall1 = wall2
  end subroutine ptock

  subroutine ptot(descr, timer, unit)
    character(len=*),     intent(in)              :: descr
    type(ptimer), target, intent(inout), optional :: timer
    integer,              intent(in),    optional :: unit

    integer               :: lun
    type(ptimer), pointer :: t
    t                      => default_timer
    if (present(timer)) t  => timer
    lun                    =  t%unit
    if (present(unit)) lun =  unit

    write(lun, t%fmt) descr, real(t%walltot, DPk)/t%countrate, &
         &                       (t% cputot)
    if (t%flush) flush(lun)
  end subroutine ptot

  subroutine ptime_print(descr, timer, unit)
    character(len=*),     intent(in)              :: descr
    type(ptimer), target, intent(inout), optional :: timer
    integer,              intent(in),    optional :: unit

    real(DPk)             :: cpu2
    integer               :: wall2, lun
    type(ptimer), pointer :: t

    t                      => default_timer
    if (present(timer)) t  => timer
    lun                    =  t%unit
    if (present(unit)) lun =  unit

    call cpu_time    (cpu2)
    call system_clock(wall2)
    t% cputot = t% cputot +  cpu2 - t% cpu1
    t%walltot = t%walltot + wall2 - t%wall1

    write(lun, t%fmt) descr, real(wall2-t%wall1, DPk)/t%countrate, &
         &                       (cpu2 -t% cpu1)
    if (t%flush) flush(lun)

    t% cpu1 =  cpu2
    t%wall1 = wall2
  end subroutine ptime_print
!!! end ptime stuff

  pure elemental function uppercase(str)
    character(*), intent(in) :: str
    character(len(str))      :: uppercase

    integer :: ic, i

    character(26), parameter :: CAP = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    character(26), parameter :: low = 'abcdefghijklmnopqrstuvwxyz'

    uppercase = str
    do i = 1, len_trim(str)
       ic = index(low, str(i:i))
       if (ic > 0) uppercase(i:i) = CAP(ic:ic)
    end do
  end function uppercase

  pure elemental function lowercase(str)
    character(*), intent(in) :: str
    character(len(str))      :: lowercase

    integer :: ic, i

    character(26), parameter :: CAP = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    character(26), parameter :: low = 'abcdefghijklmnopqrstuvwxyz'

    lowercase = str
    do i = 1, len_trim(str)
       ic = index(CAP, str(i:i))
       if (ic > 0) lowercase(i:i) = low(ic:ic)
    end do
  end function lowercase


  integer function line_count(fid)  
    !>Returns the number of lines in a file.
    !>\param fname Name of the file
    !>\return Number of lines in the file.

    implicit none

    !input parameters
    !character(len=*) fname   !filename of the file to count
    integer, intent(in) :: fid 

    !local parameters
    character(20) dummy
    integer ioStat
    logical ioEndFlag

    ioEndFlag = .false.
    line_count = 0

    do while (.not. ioEndFlag )
       read(fid,"(A20)", iostat=ioStat ) dummy
       if( iostat .eq. 0) line_count = line_count + 1
       if( iostat < 0 ) ioEndFlag = .TRUE. 
    end do
    rewind(fid)
  END FUNCTION line_count

  subroutine paropen(unit, file, iproc, iloop, iostat, status, form)
    integer,      intent(in)            :: unit
    character(*), intent(in)            :: file
    integer,      intent(in)            :: iloop
    integer,      intent(in)            :: iproc
    integer,      intent(out), optional :: iostat
    character(*), intent(in),  optional :: status
    character(*), intent(in),  optional :: form

    character(len=BUFSZ) :: tmpnam
    character(len=   11) :: f = 'FORMATTED'
    character(len=    7) :: s = 'UNKNOWN'

    if (present(status)) s = status
    if (present(form))   f = form

    if (iproc==0) then
       tmpnam = file
    else
       write(tmpnam, '(A, "_", I0)') trim(file), iloop
    end if

    if (present(iostat)) then
       open(unit, FILE=tmpnam, STATUS=s, FORM=f, IOSTAT=iostat)
    else
       open(unit, FILE=tmpnam, STATUS=s, FORM=f)
    end if
  end subroutine paropen

!!! We would like to use the NEWUNIT keyword to open(), but many
!!! compilers do not support that yet.  Workaround from
!!! <http://fortranwiki.org/fortran/show/newunit>:
!!! 
  ! This is a simple function to search for an available unit.
  ! LUN_MIN and LUN_MAX define the range of possible LUNs to check.
  ! The UNIT value is returned by the function, and also by the optional
  ! argument. This allows the function to be used directly in an OPEN
  ! statement, and optionally save the result in a local variable.
  ! If no units are available, -1 is returned.
  integer function newunit(unit)
    integer, intent(out), optional :: unit
    ! local
    integer, parameter :: LUN_MIN=1000, LUN_MAX=10000
    logical :: opened
    integer :: lun
    ! begin
    newunit=-1
    do lun=LUN_MIN,LUN_MAX
       inquire(unit=lun,opened=opened)
       if (.not. opened) then
          newunit=lun
          exit
       end if
    end do
    if (present(unit)) unit=newunit
  end function newunit
end MODULE util


!---------------  Helper procedures for command line interface --------------
module clio
  use iso_fortran_env, only: ERROR_UNIT, OUTPUT_UNIT
  use util,            only: string
  use const,           only: BUFSZ, DPk
  
  implicit none
  private

  public :: croak, carp, cluck, fetcharg, argstr

  interface fetcharg
#ifdef HAVE_VARLEN_STR
     module procedure fetcharg_alloc
#else
     module procedure fetcharg_buf
#endif

     module procedure fetcharg_str, fetcharg_int, fetcharg_real
  end interface fetcharg

  type argstr
#ifdef HAVE_VARLEN_STR
     character(len=:), allocatable :: s
#else
     character(len=BUFSZ) :: s
#endif
  end type argstr

  type(argstr) :: progname
  logical      :: got_progname = .false.
contains
  subroutine get_progname()
    integer s

    call fetcharg(0, progname, status=s)

    if (s /= 0) progname%s='WIEN2WANNIER'

    got_progname = .true.
  end subroutine get_progname

  subroutine croak(message, status)
    character(len=*), intent(in), optional   :: message
    integer,          intent(in), optional   :: status

    integer            :: s
    s=1

    if (present(status)) s=status

    if (present(message)) then
       if (.not. got_progname) call get_progname()

       write(ERROR_UNIT, '(A, ": ", A)') trim(progname%s), message
    end if

    call exit(s)
  end subroutine croak

  subroutine carp(message)
    character(len=*), intent(in) :: message

    if (.not. got_progname) call get_progname()

    write(ERROR_UNIT, '(A, ": ", A)') trim(progname%s), message
  end subroutine carp

  subroutine cluck(message)
    character(len=*), intent(in) :: message

    if (.not. got_progname) call get_progname()

    write(OUTPUT_UNIT, '(A, ": ", A)') trim(progname%s), message
  end subroutine cluck

#ifdef HAVE_VARLEN_STR
  !! This elegant variable-length version of `fetcharg´ does not work
  !! on pre-4.8 gfortran
  subroutine fetcharg_alloc(i, str, message, status)
    integer,          intent(in)            :: i
    type(argstr),     intent(out)           :: str
    character(len=*), intent(in),  optional :: message
    integer,          intent(out), optional :: status

    integer :: s, l

    if (allocated(str%s)) deallocate(str%s)

    call get_command_argument(i, length=l, status=s)
    if (present(status)) then
       status = s
    elseif (s /= 0) then
       if (present(message)) then
          call croak(message)
       else
          call croak("FETCHARG: failed to get command argument #" // &
               &     trim(string(i)) // " length: " // trim(string(s)))
       end if
    end if

    allocate(character(len=l) :: str%s)
    !! zero-length arguments seem to need special treatment
    !! (hooray Fortran!)
    if (l==0) return

    call get_command_argument(i, value=str%s, status=s)
    if (present(status)) then
       status = s
    elseif (s /= 0) then
       if (present(message)) then
          call croak(message)
       else
          call croak("FETCHARG: failed to get command argument #" // &
               &     trim(string(i)) // ": " // trim(string(s)))
       end if
    end if
  end subroutine fetcharg_alloc
#else
  subroutine fetcharg_buf(i, str, message, status)
    integer,          intent(in)               :: i
    type(argstr),     intent(out)              :: str
    character(len=*), intent(in),  optional    :: message
    integer,          intent(out), optional    :: status

    integer :: s

    !! here, a zero-length argument is okay
    !! (probably because len(str%s) /= 0)
    call get_command_argument(i, value=str%s, status=s)
    if (present(status)) then
       status = s
    elseif (s /= 0) then
       if (present(message)) then
          call croak(message)
       elseif (s < 0) then
          call croak("FETCHARG_BUF: buffer too small for command argument #"&
               &     // trim(string(i)))
       else
          call croak("FETCHARG_BUF: failed to get command argument #" // &
               &     trim(string(i)) // ": " // trim(string(s)))
       end if
    end if
  end subroutine fetcharg_buf
#endif

  subroutine fetcharg_str(i, str, message, status)
    integer,          intent(in)               :: i
    character(len=*), intent(out)              :: str
    character(len=*), intent(in),  optional    :: message
    integer,          intent(out), optional    :: status

    integer :: s

    !! here, a zero-length argument is okay
    !! (probably because len(str) /= 0)
    call get_command_argument(i, value=str, status=s)
    if (present(status)) then
       status = s
    elseif (s /= 0) then
       if (present(message)) then
          call croak(message)
       elseif (s < 0) then
          call croak("FETCHARG_STR: buffer too small for command argument #"&
               &     // trim(string(i)))
       else
          call croak("FETCHARG_STR: failed to get command argument #" // &
               &     trim(string(i)) // ": " // trim(string(s)))
       end if
    end if
  end subroutine fetcharg_str

  subroutine fetcharg_int(i, n, message, status)
    integer,          intent(in)               :: i
    integer,          intent(out)              :: n
    character(len=*), intent(in), optional     :: message
    integer,          intent(out), optional    :: status

    integer :: s
    type(argstr) :: str

    call fetcharg(i, str, STATUS=s)

    if (s /= 0) then
       if (present(status)) then
          status = s
          return
       else
          if (present(message)) then
             call croak(message)
          elseif (s < 0) then
             call croak("FETCHARG_INT: buffer too small &
                  & for command argument #" // trim(string(i)))
          else
             call croak("FETCHARG_INT: failed to get command argument #" // &
                  &     trim(string(i)) // ": " // trim(string(s)))
          end if
       end if
    end if

    read(str%s, *) n
  end subroutine fetcharg_int

  subroutine fetcharg_real(i, x, message, status)
    integer,          intent(in)               :: i
    real(DPk),        intent(out)              :: x
    character(len=*), intent(in),  optional    :: message
    integer,          intent(out), optional    :: status

    integer :: s
    type(argstr) :: str

    call fetcharg(i, str, STATUS=s)

    if (s /= 0) then
       if (present(status)) then
          status = s
          return
       else
          if (present(message)) then
             call croak(message)
          elseif (s < 0) then
             call croak("FETCHARG_INT: buffer too small &
                  & for command argument #" // trim(string(i)))
          else
             call croak("FETCHARG_INT: failed to get command argument #" // &
                  &     trim(string(i)) // ": " // trim(string(s)))
          end if
       end if
    end if

    read(str%s, *) x
  end subroutine fetcharg_real
end module clio


!---------------- ‘struct_t’ type and associated procedures   ---------------
!
! NB: individual wien2wannier programs may still use their own
! struct-related modules
module structmod
  use const, only: DPk

  implicit none

  private
  public :: struct_t, struct_read

  interface struct_read
     module procedure struct_read_unit, struct_read_fname
  end interface struct_read

  type intptr
     integer, allocatable :: p(:)
  end type intptr

  type struct_t
     !! Structure to represent a ‘struct’ file
     !!
     !! The conventions here follow mostly those of ‘structeditor’

     character(len=80) :: title
     character(len= 4) :: lattic       ! lattice type
     integer           :: nneq, nat    ! number of noneq. and total atoms
     character(len= 4) :: mode         ! ‘RELA’ or ‘NREL’
     real(DPk)         :: a(3)         ! lattice constants
     real(DPk)         :: alpha(3)     ! angles
     real(DPk)         :: gbas(3,3)    ! reciprocal lattice vectors
     real(DPk)         :: brlat(3,3)   ! Bravais lattice (row is a vector)
     real(DPk)         :: lat2car(3,3) ! lattice to cartesian transformation
                                       ! (different from brlat in some cases)
     real(DPk)         :: stru2frac(3,3)!‘struct’ to “fractional” coords. 
     logical           :: ortho, prim  ! orthogonal axes?  primitive latt.?
     real(DPk)         :: vol          ! u.c. volume
     integer, allocatable, dimension(:,:,:) &
          &            :: rsym, ksym   ! symmetry operations (3×3×N — can be
                                       ! different!) as ‘outputkgen’
     real(DPk), allocatable, dimension(:,:) &
          &            :: rtrans       ! translational part of symmetries
                                       ! (3×N)

     ! positions (3 × nat); local rotation matrices (3 × 3 × nneq)
     real(DPk),         allocatable :: pos(:,:), locrot(:,:,:)
     ! neq2at(ineq) is the atom numer corresponding to ineq
     integer,           allocatable :: mult(:), isplit(:), npt(:), neq2at(:)
     ! at2neq(iat)%p is an array of the ‘ineq’s corresponding to iat
     type(intptr),      allocatable :: at2neq(:)
     character(len=10), allocatable :: aname(:)
     real(DPk),         allocatable :: r0(:), rmt(:), Z(:)
  end type struct_t
contains
  subroutine struct_read_unit(lun, stru)
    use const, only: PI, SQ3, ORTHO_TEST
    use clio,  only: croak
    use util,  only: inverse3x3

    integer,        intent(in)  :: lun
    type(struct_t), intent(out) :: stru

    integer   :: iat, ineq, N
    ! “Bravais lattice” and “lattice->cartesian” transformation
    ! matrices [in some cases, the two are different], both for direct
    ! and reciprocal space
    real(DPk) :: br1_dir(3,3), br1_rec(3,3), br2_dir(3,3), br2_rec(3,3)
    ! abbreviations
    real(DPk) :: pia(3), alpha(3), cosab, cosac, cosbc, sinab, sinbc
    real(DPk) :: rvfac, wurzel

    read(lun, '(A)')          stru%title
    read(lun, '(A, 23X, I3)') stru%lattic, stru%nneq
    read(lun, '(13X, A)')     stru%mode
    read(lun, '(6F10.6)')     stru%a, stru%alpha

    where (abs(stru%alpha) < 1e-5) stru%alpha = 90

    N = stru%nneq
    allocate(stru%mult(N), stru%isplit(N), stru%aname (N), &
         &   stru%npt (N), stru%at2neq(N), stru%Z     (N), &
         &   stru%r0  (N), stru%rmt   (N), stru%locrot(3,3,N))

    countat: do iat = 1, stru%nneq
       read(lun,*) ! first position
       read(lun, '(15X, I2, 17X, I2)') stru%mult(iat), stru%isplit(iat)
       ! rest of positions
       do ineq = 2, stru%mult(iat)
          read(lun,*)
       end do
       read(lun, '(A, 5X, I5, 5X, F10.8, 5X, F10.5, 5X, F5.2)') &
            stru%aname(iat), stru%npt(iat), stru%r0(iat), stru%rmt(iat), &
            stru%Z(iat)
       read(lun,'(20X, 3F10.7)') stru%locrot(1, :, iat)
       read(lun,'(20X, 3F10.7)') stru%locrot(2, :, iat)
       read(lun,'(20X, 3F10.7)') stru%locrot(3, :, iat)
    end do countat

    stru%nat = sum(stru%mult)
    N = stru%nat
    allocate(stru%pos(3, N), stru%neq2at(N))

    N=0
    do iat = 1, stru%nneq
       stru%neq2at(N+1:N+stru%mult(iat)) = iat

       allocate(stru%at2neq(iat)%p(stru%mult(iat)))
       stru%at2neq(iat)%p = iat

       N = N+stru%mult(iat)
    end do

    rewind(lun)
    read(lun,*) ! title
    read(lun,*) ! lattic
    read(lun,*) ! mode
    read(lun,*) ! uc
    N=0
    readpos: do iat = 1, stru%nneq
       N=N+1
       read(lun, '(9X, 3(3X, F10.8))') stru%pos(:, N)
       read(lun,*) ! mult
       do ineq = 2, stru%mult(iat)
          N=N+1
          read(lun, '(9X, 3(3X, F10.8))') stru%pos(:, N)
       end do
       read(lun,*) ! aname
       read(lun,*) ! locrot
       read(lun,*) !
       read(lun,*) !
    end do readpos

    !! Now compute lattice vectors.  This is copied from
    !! SRC_structeditor/module.f (subroutine latgen_struct)
    pia = 2*PI/stru%a
    alpha = stru%alpha*PI/180

    cosab = cos(alpha(3)); sinab = sin(alpha(3))
    cosac = cos(alpha(2))
    cosbc = cos(alpha(1)); sinbc = sin(alpha(1))

    br1_rec = 0; br1_dir = 0; br2_rec = 0; br2_dir = 0

    lattic: select case (stru%lattic(1:1))
    case ('H')
       br1_rec(1,1) = pia(1)*2/SQ3
       br1_rec(1,2) = pia(1)  /SQ3
       br1_rec(2,2) = pia(2)
       br1_rec(3,3) = pia(3)

       br2_rec = br1_rec

       rvfac = 2/SQ3

       stru%ortho = .false.
       stru%prim  = .true.

    case ('S', 'P')             ! what's ‘S’??
       wurzel = sqrt(sinbc**2 - cosac**2 - cosab**2 + 2*cosbc*cosac*cosab)

       br1_rec(1,1) = pia(1) * sinbc/wurzel
       br1_rec(1,2) = pia(2) * (-cosab + cosbc*cosac)/(sinbc*wurzel)
       br1_rec(1,3) = pia(3) * (-cosac + cosbc*cosab)/(sinbc*wurzel)
       br1_rec(2,2) = pia(2) / sinbc
       br1_rec(2,3) =-pia(3) * cosbc/sinbc
       br1_rec(3,3) = pia(3)

       br2_rec = br1_rec

       rvfac = 1/wurzel

       stru%ortho = all(abs(alpha - PI/2) <= ortho_test)
       stru%prim  = .true.

    case ('F')
       br1_rec(1,1) = pia(1)
       br1_rec(2,2) = pia(2)
       br1_rec(3,3) = pia(3)

       br2_rec(1,:) = pia(1) * (/-1, 1, 1 /)
       br2_rec(2,:) = pia(2) * (/ 1,-1, 1 /)
       br2_rec(3,:) = pia(3) * (/ 1, 1,-1 /)

       rvfac = 4
       stru%ortho = .true.
       stru%prim  = .false.

    case ('B')
       br1_rec(1,1) = pia(1)
       br1_rec(2,2) = pia(2)
       br1_rec(3,3) = pia(3)

       br2_rec(1,:) = pia(1) * (/ 0, 1, 1 /)
       br2_rec(2,:) = pia(2) * (/ 1, 0, 1 /)
       br2_rec(3,:) = pia(3) * (/ 1, 1, 0 /)

       rvfac = 2

       stru%ortho = .true.
       stru%prim  = .false.

    case ('R')
       br1_rec(1, :) = pia(1) * (/  1, 1, -2 /)/SQ3
       br1_rec(2, :) = pia(2) * (/ -1, 1,  0 /)
       br1_rec(3, :) = pia(3)

       br2_rec = br1_rec

       rvfac = 6/SQ3

       stru%ortho = .false.
       stru%prim  = .true.

    case ('C')
       ! “defaults”, to be changed in nonorthogonal XZ-case
       br1_rec(1,1) = pia(1)
       br1_rec(2,2) = pia(2)
       br1_rec(3,3) = pia(3)

       rvfac = 2
       stru%ortho = .true.      ! may be overridden for CXZ
       stru%prim  = .false.

       C: select case (stru%lattic(2:3))
       case ('XY')
          br2_rec(1,:) = pia(1) * (/ 1, 1, 0 /)
          br2_rec(2,:) = pia(2) * (/-1, 1, 0 /)
          br2_rec(3,:) = pia(3) * (/ 0, 0, 1 /)

       case ('YZ')
          br2_rec(1,:) = pia(1) * (/ 1, 0, 0 /)
          br2_rec(2,:) = pia(2) * (/ 0, 1, 1 /)
          br2_rec(3,:) = pia(3) * (/ 0,-1, 1 /)

       case ('XZ')
          ortho: if (abs(alpha(3) - PI/2) <= 0.0001) then
             br2_rec(1,:) = pia(1) * (/ 1, 0, 1 /)
             br2_rec(2,:) = pia(2) * (/ 0, 1, 0 /)
             br2_rec(3,:) = pia(3) * (/-1, 0, 1 /)
          else                  ! CXZ monoclinic case
             br1_rec(1,1) = pia(1) / sinab
             br1_rec(1,2) =-pia(2) / sinab * cosab

             br2_rec(1,1) = pia(1) / sinab
             br2_rec(1,2) =-pia(2) / sinab * cosab
             br2_rec(1,3) = pia(1) / sinab
             br2_rec(2,:) = pia(2) * (/ 0, 1, 0 /)
             br2_rec(3,:) = pia(3) * (/-1, 0, 1 /)

             rvfac = 2/sinab
             stru%ortho = .false.
          end if ortho
       end select C

    case default
       rvfac = 0                ! silence warning
       call croak('unknown lattice type `' // trim(stru%lattic) // "'")
    end select lattic

    call inverse3x3(br1_rec, br1_dir);
    call inverse3x3(br2_rec, br2_dir);

    stru%vol = product(stru%a) / rvfac
    stru%gbas      = br2_rec
    stru%lat2car   = br1_dir*2*PI
    stru%brlat     = br2_dir*2*PI
    stru%stru2frac = matmul(br1_dir, br2_rec)

    ! make_ksym() needs gbas and friends
    call read_rsym()
    call make_ksym()
  contains
    subroutine read_rsym()
!!! This is from SRC_kgen (main.f, addinv.f)
      integer :: Nrsym, i, id

      read(lun, '(I4)') Nrsym
      allocate(stru%rsym  (3,3,Nrsym), &
           &   stru%rtrans(3,  Nrsym))
      do i=1,Nrsym
         read(lun, '(3(3I2,F11.8/), I8)') &
              (stru%rsym(:,id,i), stru%rtrans(id, i), id=1,3)
      end do
    end subroutine read_rsym

    subroutine make_ksym()
!!! This is from SRC_kgen/sdefl.f

      real(DPk), dimension(3,3) :: a
      logical                   :: addinv
      integer                   :: i, Nrsym, Nksym

      addinv = .true.
      Nrsym  = size(stru%rsym, 3)
!!! Redefine symmetry operations from struct with
!!! unitary transformation  u(-1) * S * u

      do i=1,Nrsym
         ! if inversion is already present, we do not need to add it
         ! later
         if (all(stru%rsym(:,:,i) == reshape( (/ -1, 0, 0,    &
              &                                   0,-1, 0,    &
              &                                   0, 0,-1 /), &
              &                               (/ 3, 3 /) ))) then
            addinv = .false.
            exit
         end if
      end do

      if (Nrsym>24 .and. addinv) &
           call croak('Nsym > 24 without inversion')

      if (addinv) then
         ! add inversion
         ! FIXME: check if we should actually do this (non-sp non-so)
         Nksym = 2*Nrsym
      else
         Nksym =   Nrsym
      end if
      allocate(stru%ksym(3,3,Nksym))

                  stru%ksym(:,:,       1:Nrsym) =  stru%rsym(:,:, 1:Nrsym)
      if (addinv) stru%ksym(:,:, Nrsym+1:Nksym) = -stru%rsym(:,:, 1:Nrsym)

      if (stru%ortho .or. stru%lattic=='CXZ') then
         ! FIXME: why the check for CXZ??  sdefl.f has it …
         do i=1,Nksym
            a = matmul(matmul(stru%brlat, stru%ksym(:,:,i)), &
                 &     stru%gbas)/2/PI
            stru%ksym(:,:,i) = nint(a)
         end do
      end if
    end subroutine make_ksym
  end subroutine struct_read_unit

  subroutine struct_read_fname(fname, stru)
    use util, only: newunit

    character(*),   intent(in)  :: fname
    type(struct_t), intent(out) :: stru

    integer :: lun

    open(unit=newunit(lun), file=fname, status='old')
    call struct_read_unit(lun, stru)
    close(lun)
  end subroutine struct_read_fname
end module structmod


!---------------- ‘inwf_t’ type and associated procedures   -----------------
module inwfmod
  use const, only: DPk

  implicit none

  private
  public :: inwf_t, inwf_read

  interface inwf_read
     module procedure inwf_read_unit, inwf_read_fname, inwf_read_argstr
  end interface inwf_read

  type proj
     integer                                 :: NY
     integer,      allocatable, dimension(:) :: iat, l, m
     complex(DPk), allocatable, dimension(:) :: coeff
  end type proj

  type inwf_t
     !! Structure to represent a ‘inwf’ file
     character(len=4)    :: mode
     logical             :: Mmn, Amn
     integer             :: bmin, bmax, Nproj, LJmax
     type(proj), pointer :: projections(:)
  end type inwf_t

contains
  subroutine inwf_read_unit(lun, inwf)
    use clio,  only: croak
    use util,  only: uppercase

    integer,      intent(in)  :: lun
    type(inwf_t), intent(out) :: inwf

    integer             :: w, y, NY
    real(DPk)           :: re, im
    type(proj), pointer :: p

    read(lun, *) inwf%mode

    select case (uppercase(inwf%mode))
    case('AMN')
       inwf%amn = .true.
    case('MMN')
       inwf%mmn = .true.
    case('BOTH')
       inwf%amn = .true.
       inwf%mmn = .true.
    case default
       call croak('unknown MODE in inwf, '//inwf%mode)
    end select
       
    read(lun, *) inwf%bmin, inwf%bmax
    if (inwf%bmax < inwf%bmin) &
         call croak('bmin > bmax in inwf')

    read(lun,*) inwf%LJmax, inwf%Nproj

    if (inwf%Amn) then
       allocate(inwf%projections(inwf%Nproj))
       do w = 1, inwf%Nproj
          read(lun, *) NY
          p => inwf%projections(w)
          p%NY = NY
          allocate(p%iat(NY), p%l(NY), p%m(NY), p%coeff(NY))

          ylm: do y = 1,NY
             read(lun, *) p%iat(y), p%l(y), p%m(y), re, im
             p%coeff(y) = cmplx(re, im, DPk)
          end do ylm
       end do
    end if
  end subroutine inwf_read_unit

  subroutine inwf_read_fname(fname, inwf)
    use util, only: newunit

    character(*), intent(in)  :: fname
    type(inwf_t), intent(out) :: inwf

    integer :: lun

    open(unit=newunit(lun), file=fname, status='old')
    call inwf_read_unit(lun, inwf)
    close(lun)
  end subroutine inwf_read_fname

  subroutine inwf_read_argstr(arg, inwf)
    use clio, only: argstr

    type(argstr), intent(in)  :: arg
    type(inwf_t), intent(out) :: inwf

    call inwf_read_fname(arg%s, inwf)
  end subroutine inwf_read_argstr
end module inwfmod


!---------------- Reading ‘klist’ files            --------------------------
module kpoints
  use const,     only: BUFSZ, DPk
  use structmod, only: struct_t
  use util,      only: string, newunit
  use clio,      only: croak

  implicit none
  private
  public :: get_kmesh_klist, get_kmesh_band, &
       &    count_kmesh_klist, count_kmesh_band

  interface get_kmesh_klist
     module procedure get_kmesh_klist_unit, get_kmesh_klist_fname
  end interface get_kmesh_klist

  interface get_kmesh_band
     module procedure get_kmesh_band_unit, get_kmesh_band_fname
  end interface get_kmesh_band

  interface count_kmesh_klist
     module procedure count_kmesh_klist_unit, count_kmesh_klist_fname
  end interface count_kmesh_klist

  interface count_kmesh_band
     module procedure count_kmesh_band_unit, count_kmesh_band_fname
  end interface count_kmesh_band

  character(len=*), parameter :: &
       fmt_klist_old   = "(A10,4I 5)",         &
       fmt_klist_new   = "(A10,4I10)",         &
       fmt_klist_old_1 = "(10X,4I 5,35X,3I3)", &
       fmt_klist_new_1 = "(10X,4I10,35X,3I3)"
  integer,      parameter :: fmt_klist_len = len(fmt_klist_new)
contains
!!! Read standard ‘klist’
  subroutine get_kmesh_klist_unit(unit, kpoints, stru, mp_grid)
    integer,        intent(in)               :: unit
    real(DPk),      intent(out), allocatable :: kpoints(:,:)
    type(struct_t), intent(in)               :: stru
    integer,        intent(out), optional    :: mp_grid(3)

    integer                  :: nk, ik
    character(fmt_klist_len) :: fmt
    character(10)            :: kname
    integer                  :: kfrac(3), div
    real(DPk)                :: RR(3,3)

    if (stru%ortho) then
       RR(1,:) = stru%brlat(1,:) / stru%a
       RR(2,:) = stru%brlat(2,:) / stru%a
       RR(3,:) = stru%brlat(3,:) / stru%a
    end if

    if (present(mp_grid)) then
       call klist_detect_fmt(unit, fmt, mp_grid)
    else
       call klist_detect_fmt(unit, fmt)
    end if

    nk = count_kmesh_klist(unit)
    allocate(kpoints(nk, 3))

    ik = 0
    readk: do
       read(unit, fmt) kname, kfrac, div
       if (kname == "END") exit

       ik = ik+1

       kpoints(ik,:) = real(kfrac, DPk)/real(div, DPk)

       if (stru%ortho) kpoints(ik,:) = matmul(RR, kpoints(ik,:))
    end do readk

    if (ik /= nk) &
         call croak('get_kmesh_klist: numbers of k-points do not match' // &
         &          trim(string(ik)) // ' /= ' // trim(string(nk)))
  end subroutine get_kmesh_klist_unit

!!! Read ‘klist_band’, taking into account only labeled points
  subroutine get_kmesh_band_unit(unit, kpoints, stru, knames)
!!! Reading the band-structure klist gets its own subroutine because
!!! the required work (k-point labels, #points per segment) as well as
!!! the structure of the file (no running numbers, no mp_grid, ...)
!!! are sufficiently different.

    integer,                intent(in)  :: unit
    real(DPk), allocatable, intent(out) :: kpoints(:,:)
    type(struct_t),         intent(in)  :: stru
    character, allocatable, intent(out) :: knames(:)

    integer                  :: nk, ik
    character(fmt_klist_len) :: fmt
    character(10)            :: kname
    integer                  :: kfrac(3), div
    real(DPk)                :: RR(3,3)

    if (stru%ortho) then
       RR(1,:) = stru%brlat(1,:) / stru%a
       RR(2,:) = stru%brlat(2,:) / stru%a
       RR(3,:) = stru%brlat(3,:) / stru%a
    end if

    call klist_detect_fmt(unit, fmt)
    nk = count_kmesh_band(unit)

    allocate(kpoints(nk, 3), knames(nk))

    ik  = 0
    readk: do
       read(unit, fmt) kname, kfrac, div
       if (kname == "END") exit
       if (kname == "")    cycle

       ! now we have a k-vector we should save
       ik = ik+1
       knames(ik) = kname(1:1)

       kpoints(ik,:) = real(kfrac, DPk)/real(div, DPk)

       if (stru%ortho) &
            kpoints(ik,:) = matmul(RR, kpoints(ik,:))
    end do readk

    if (ik /= nk) &
         call croak('get_kmesh_band: numbers of k-points do not match' // &
         &          trim(string(ik)) // ' /= ' // trim(string(nk)))
  end subroutine get_kmesh_band_unit

!!! Count k-points in standard ‘klist’
  function count_kmesh_klist_unit(unit) result(nk)
    integer, intent(in)  :: unit
    integer              :: nk

    character(len=BUFSZ) :: buf

    nk = 0
    countk: do
       read(unit, '(A)') buf
       if (buf == "END") exit
       nk = nk+1
    end do countk
    rewind(unit)
  end function count_kmesh_klist_unit

!!! Count labelled k-points in ‘klist_band’
  function count_kmesh_band_unit(unit) result(nk)
    integer, intent(in)  :: unit
    integer              :: nk

    character(len=6)     :: buf

    nk = 0
    countk: do
       read(unit, '(A)') buf
       if (buf == "END") exit
       ! we only count k-points with a label, assuming that these
       ! delimit "paths" that will be passed to wannier90
       if (buf /= "") nk = nk+1
    end do countk
    rewind(unit)
  end function count_kmesh_band_unit

!!! Wrappers for calling the above with a filename instead of an open unit
  subroutine get_kmesh_klist_fname(fname, kpoints, stru, mp_grid)
    character(*),   intent(in)               :: fname
    real(DPk),      intent(out), allocatable :: kpoints(:,:)
    type(struct_t), intent(in)               :: stru
    integer,        intent(out), optional    :: mp_grid(3)

    integer :: lun

    open(unit=newunit(lun), file=fname, status='old')
    call get_kmesh_klist_unit(lun, kpoints, stru, mp_grid)
    close(lun)
  end subroutine get_kmesh_klist_fname

  subroutine get_kmesh_band_fname(fname, kpoints, stru, knames)
    character(*),           intent(in)  :: fname
    real(DPk), allocatable, intent(out) :: kpoints(:,:)
    type(struct_t),         intent(in)  :: stru
    character, allocatable, intent(out) :: knames(:)

    integer :: lun

    open(unit=newunit(lun), file=fname, status='old')
    call get_kmesh_band_unit(lun, kpoints, stru, knames)
    close(lun)
  end subroutine get_kmesh_band_fname

  function count_kmesh_klist_fname(fname) result(nk)
    character(*), intent(in)  :: fname
    integer                   :: nk

    integer :: lun

    open(unit=newunit(lun), file=fname, status='old')
    nk = count_kmesh_klist_unit(lun)
    close(lun)
  end function count_kmesh_klist_fname

  function count_kmesh_band_fname(fname) result(nk)
    character(*), intent(in)  :: fname
    integer                   :: nk

    integer :: lun

    open(unit=newunit(lun), file=fname, status='old')
    nk = count_kmesh_band_unit(lun)
    close(lun)
  end function count_kmesh_band_fname

!!! Try to figure out which format our klist file follows, optionally
!!! reading ‘mp_grid’ in the process.
  subroutine klist_detect_fmt(unit, fmt, mp_grid)
    integer,                  intent(in)  :: unit
    character(fmt_klist_len), intent(out) :: fmt
    integer,        optional, intent(out) :: mp_grid(3)

    character(BUFSZ) :: iom
    integer          :: idum(4), ios

    ! First, we try the new format:
    if (present(mp_grid)) then
       read(unit, fmt_klist_new_1, iostat=ios, iomsg=iom) idum, mp_grid
    else
       read(unit, fmt_klist_new_1, iostat=ios, iomsg=iom) idum
    end if
    if (ios == 0) then
       ! Success — new format it is
       fmt = fmt_klist_new
    else if (ios < 0) then
       call croak('get_kmesh_klist(): EOF encountered while reading first&
            & line of ‘klist’ file: '//trim(iom))
    else
       ! Error in read, try old format
       if (present(mp_grid)) then
          read(unit, fmt_klist_old_1, iostat=ios, iomsg=iom) idum, mp_grid
       else
          read(unit, fmt_klist_old_1, iostat=ios, iomsg=iom) idum
       end if
       if (ios == 0) then
          fmt = fmt_klist_old
       else
          call croak('get_kmesh_klist(): something went wrong while reading&
               & first line of ‘klist’ file: '//trim(iom))
       end if
    end if

    rewind(unit)
  end subroutine klist_detect_fmt
end module kpoints


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


!--------------- Procedures nicked from Wien2k       ------------------------
module wien2k
  private
  public :: errflg, errclr, gtfnam
  contains
      SUBROUTINE ERRFLG(FNAME,MSG)
      CHARACTER(*)      FNAME, MSG
!
!     ..................................................................
!
! 1.     PROGRAM UNIT 'ERRFLG'
!           Notify that an error has occured.
!           FORTRAN 77 SUBROUTINE
!
! 2.     PURPOSE
!           Because there is no standard (or even semi-standard) way to
!           generate exit codes in FORTRAN 77, this routine writes a
!           non-empty file to the current subdirectory as an indication
!           that some serious error has occured. Other programs can then
!           check the contents of this file to determine whether an
!           error has occured. The errorflag-file is left opened when
!           returning from this routine to enable writing other
!           errormessages to it.
!
! 3.     USAGE
!           CALL ERRFLG('lapw2.error','Error in OUTWIN')
!
!        ARGUMENT-DESCRIPTION
!           FNAME  - CHARACTER*(*) string                        (input)
!                    The name of the file acting as error-flag.
!           MSG    - CHARACTER*(*) string                        (input)
!                    The (error) message which should be written to the
!                    errorflag-file.
!
!        USED SUBROUTINES (DIRECTLY CALLED)
!           none
!
!        INDIRECTLY CALLED SUBROUTINES
!           none
!
!        UTILITY-SUBROUTINES (USE BEFOREHAND OR AFTERWARDS)
!           ERRCLR - clears the contents of a file
!
!        INPUT/OUTPUT (READ/WRITE)
!           A message given as argument MSG is written to a the file
!           given as argument FNAME. File FNAME is created if not
!           existing and otherwise overwritten.
!
!        MACHINENDEPENDENT PROGRAMPARTS
!           none
!
! 4.     REMARKS
!           The best way to use this routine is to call ERRFLG at the
!           start of a program writing some message to the errorflag-
!           file and ERRCLR before a successful exit of the program. By
!           checking the contents of the errorflag-file it is possible
!           to determine if the program was successfully completed.
!
!           This method has the advantage of working even if some
!           runtime-error occurs which is not taken care of in the
!           program.
!
! 5.     METHOD
!           - open errorflag-file
!           - write some message to this errorflag-file
!           - exit leaving the file opened
!
! 6.     DATE
!           24. August 1993                                 Version 1.02
!
!        INSTITUT FUER TECHNISCHE ELEKTROCHEMIE            --  TU VIENNA
!     ..................................................................
!
      OPEN (99,FILE=FNAME,ERR=900)
      WRITE (99,9000) MSG
      CLOSE (99)
      OPEN (99,FILE=FNAME,ERR=900)
!
      RETURN
!
!        Errors
!
  900 STOP 'ERRFLG - couldn''t open errorflag-file.'
!
 9000 FORMAT (A)
!
!        End of 'ERRFLG'
!
    END SUBROUTINE ERRFLG

    
    SUBROUTINE ERRCLR(FNAME)
      CHARACTER*(*)      FNAME

      CLOSE (99)
      OPEN (99,FILE=FNAME,STATUS='REPLACE')
      RETURN
    END SUBROUTINE ERRCLR

      SUBROUTINE GTFNAM(DEFFN,ERRFN,IPROC)
      implicit integer (a-n)
      CHARACTER(*)      DEFFN, ERRFN
      INTEGER            IPROC
!
!     ..................................................................
!
! 1.     PROGRAM UNIT 'GTFNAM'
!           Provide names for definition- and error-files for LAPW.
!           FORTRAN 77 SUBROUTINE
!
! 2.     PURPOSE
!           Read the commandline-argument 
!           specifying the name of the definition file ('lapw1.def')
!           and generate the name of the error file by replacing the
!           extension of the definition filename with '.error'. If no
!           extension can be found '.error' is appended.
!
!           For the parallel version in the second commandline parameter
!           the number of parallel processes is specified
!
! 3.     USAGE
!           CHARACTER*80 DEFFN, ERRFN
!           CALL GTFNAM(DEFFN,ERRFN)
!
!        ARGUMENT-DESCRIPTION
!           DEFFN  - CHARACTER*(*) string                       (output)
!                    on exit contains the filename of the lapw2-
!                    definition file 'lapw2.def'.
!           ERRFN  - CHARACTER*(*) string                       (output)
!                    on exit contains the filename of the file where
!                    error messages are stored (derived from the file-
!                    name of the definition-file).
!           IPROC  - number of parallel processes (if specified)
!
!
!        USED SUBROUTINES (DIRECTLY CALLED)
!           none
!
!        INDIRECTLY CALLED SUBROUTINES
!           none
!
!        UTILITY-SUBROUTINES (USE BEFOREHAND OR AFTERWARDS)
!           none
!
!        INPUT/OUTPUT (READ/WRITE)
!           none
!
!        MACHINENDEPENDENT PROGRAMPARTS
!           - Subroutine GETARG for the extraction of command-line
!             arguments is used (the index for referencing command-
!             line argument starts with 0 on some machines and with
!             1 on other machines).
!           - A compiler directive to enable the use of the extension
!             subroutine GETARG is used (HP-Version)
!
! 4.     REMARKS
!           It is assumed that filename-extensions are separated by
!           character '.'.
!
! 5.     METHOD
!           - get commandline-argument (taking into account that
!             comandline-arguments are referenced starting with index
!             0 on some machines and starting with index 1 on others).
!           - 'lapw1.def' := commandline-argument
!           - look for the last occurence of character '.' in the
!             commandline-argument
!           - if found replace all characters after that '.' with
!             'error' giving the error filename
!           - otherwise append '.error' giving the error filename
!           
! 6.     DATE
!           25. August 1993                                 Version 1.01
!
!        INSTITUT FUER TECHNISCHE ELEKTROCHEMIE            --  TU VIENNA
!     ..................................................................
!
!        Local Parameters
      !
      CHARACTER(5)       ERREXT
      PARAMETER          (ERREXT = 'error')
!
!        Local Scalars
!
      INTEGER            I, iarg
!
!        extract the command-line argument
!
      IPROC=0
      iarg=command_argument_count()
      if(iarg.eq.1) then
         CALL get_command_argument(iarg,DEFFN)
      else if(iarg.eq.2) then
         CALL get_command_argument(2,DEFFN)
         READ(DEFFN,*)IPROC
         CALL get_command_argument(1,DEFFN)
      else
         STOP 'GTFNAM - One or two commandline arguments have to be given.'
      endif
!
!        generate a name for the error-message file
!
      DO 10 I = LEN(DEFFN), 1, -1
         IF (DEFFN(I:I) .EQ. '.') THEN
            IF (LEN(ERRFN) .LT. (I+LEN(ERREXT))) GOTO 910
            ERRFN(1:I) = DEFFN(1:I)
            ERRFN(I+1:LEN(ERRFN)) = ERREXT
            GOTO 30
         ENDIF
   10 CONTINUE
!
!        the name of the definition file contains no '.', it is assumed
!        that this name contains no extension - append the extension
!        '.error' to get a name for the error file.
!
      DO 20 I = LEN(DEFFN), 1, -1
         IF (DEFFN(I:I) .NE. ' ') THEN
            IF (LEN(ERRFN) .LT. (I+1+LEN(ERREXT))) GOTO 910
            ERRFN(1:I) = DEFFN(1:I)
            ERRFN(I+1:LEN(ERRFN)) = '.' // ERREXT
            GOTO 30
         ENDIF
   20 CONTINUE
!
!        filename contains only spaces
!
      STOP 'GTFNAM - string ERRFN contains just spaces.'
   30 CONTINUE
!
      RETURN
!
!        Errors
!
  910 STOP 'GTFNAM - string ERRFN too short to hold filename.'
!
!        End of 'GTFNAM'
!
      END SUBROUTINE
END MODULE wien2k


!!/---
!! Local Variables:
!! mode: f90
!! End:
!!\---
!!
!! Time-stamp: <2015-12-23 12:34:19 assman@faepop36.tu-graz.ac.at>
