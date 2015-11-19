!!! wien2wannier/SRC_w2w/util_w2w.F
!!!
!!!    Collection of routines for wien2wannier & woptic
!!!
!!! Copyright 2010-2012 Philipp Wissgott
!!!           2013-2014 Elias Assmann
!!!
!!! $Id: util_w2w.F 285 2014-10-06 16:14:51Z assmann $


!!/=== Things that are exported here =============
!!
!! const: Constants mathematical and configurational
!!
!!    BUFSZ, DPk, IPk, R4, R8, C16
!!    PI, ORTHO_TEST, SQ3
!!
!!
!! reallocate: Reallocation of various arrays
!!
!!    realloc()
!!
!!
!! util: Assorted miscellanea
!!
!!    string(), uppercase(), lowercase(),
!!    inverse3x3(), ptime(), line_count(), paropen(), newunit()
!!
!!
!! clio: Helper procedures for command line interface
!!
!!    croak(), carp(), cluck(), fetcharg(), type(argstr)
!!
!!
!! structmod: ‘struct’ type and associated procedures
!!
!!    type(struct), struct_read()
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
!!    type(chk_file), chk_read()
!!
!!
!! woptic: Things specific to woptic
!!
!!    get_mommat(), get_mommat_1k(), type(inwop), inwop_read()
!!
!!
!! woptic_io: Files and units for woptic
!!
!!    fn_*, unit_*, suf_*; print1or3(), set_casename()
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

  integer, parameter  :: BUFSZ = 256
  !! Default kinds
  integer, parameter  :: DPk = kind(1.0d0) ! inherited from W90 2.0.0
  integer, parameter  :: IPk = selected_int_kind(8)
  !! Kinds for WIEN-compatibility
  real*4,     private :: four_real
  real*8,     private :: eight_real
  complex*16, private :: sixteen_complex
  integer, parameter  :: R4 = kind(four_real)
  integer, parameter  :: R8 = kind(eight_real)
  integer, parameter  :: C16= kind(sixteen_complex)

  real(DPk), parameter:: PI         = 3.1415926535897932_DPk
  real(DPk), parameter:: ORTHO_TEST = 1.e-12_DPk
  real(DPk), parameter:: SQ3        = sqrt(3._DPk)
  real(DPk), parameter:: Ryd_eV     = 13.605698_DPk ! value for WIEN compat.
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
  use const, only: DPk, BUFSZ
  implicit none

  private
  public :: string, inverse3x3, ptime, uppercase, lowercase, line_count, &
       &    paropen, newunit

  interface string
     module procedure int2str, real2str
  end interface string
contains
  character(len=10) function int2str(n)
    integer, intent(in) :: n
    write(int2str,"(I0)") n
  end function int2str

  character(len=15) function real2str(r)
    real(DPk), intent(in) :: r
    write(real2str,"(E15.8)") r
  end function real2str

  subroutine inverse3x3(a, ainv)
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

  subroutine ptime(descr, unit)
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
    use iso_fortran_env, only: ERROR_UNIT
    implicit none

    character(len=*), intent(in), optional :: descr
    integer,          intent(in), optional :: unit

    character(len=*), parameter :: fmt = "('Times for ', A, T33, '(sec):', &
         & F8.3, ' wall;', F9.3)"!, ' cpu =', F8.3, ' user +', F8.3, ' sys')"

    real(DPk),  save :: cputime1, cputime2
    integer,    save :: walltime1, walltime2, count_rate
    integer,    save :: default_lun=ERROR_UNIT
    integer          :: lun

    if (.not. present(descr)) then
       call cpu_time(cputime1)
       call system_clock(walltime1, count_rate)

       if (present(unit)) default_lun=unit

       return
    end if

    if (present(unit)) then
       lun=unit
    else
       lun=default_lun
    end if

    call cpu_time(cputime2)
    call system_clock(walltime2)

    write(lun, fmt) descr, real(walltime2-walltime1)/count_rate, &
         & (cputime2-cputime1)

    walltime1 = walltime2
    cputime1  = cputime2
  end subroutine ptime

  pure function uppercase(str)
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

  pure function lowercase(str)
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
    integer, parameter :: LUN_MIN=10, LUN_MAX=1000
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
  !! on pre-4.8 gcc
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


!---------------- ‘struct’ type and associated procedures   -----------------
!
! NB: individual wien2wannier programs may still use their own
! struct-related modules
module structmod
  use const, only: DPk

  implicit none

  private
  public :: struct, struct_read

  interface struct_read
     module procedure struct_read_unit, struct_read_fname
  end interface struct_read

  type intptr
     integer, allocatable :: p(:)
  end type intptr

  type struct
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
     logical           :: ortho
     real(DPk)         :: vol          ! u.c. volume
     integer, allocatable, dimension(:,:,:) :: rsym, ksym  ! symmetry operations from struct
                                       ! (as ‘outputkgen’)

     ! positions (3 × nat); local rotation matrices (3 × 3 × nneq)
     real(DPk),         allocatable :: pos(:,:), locrot(:,:,:)
     ! neq2at(ineq) is the atom numer corresponding to ineq
     integer,           allocatable :: mult(:), isplit(:), npt(:), neq2at(:)
     ! at2neq(iat)%p is an array of the ‘ineq’s corresponding to iat
     type(intptr),      allocatable :: at2neq(:)
     character(len=10), allocatable :: aname(:)
     real(DPk),         allocatable :: r0(:), rmt(:), Z(:)
  end type struct
contains
  subroutine struct_read_unit(lun, stru)
    use const, only: PI, SQ3, ORTHO_TEST
    use clio,  only: croak
    use util,  only: inverse3x3

    integer,      intent(in)  :: lun
    type(struct), intent(out) :: stru

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

    case ('F')
       br1_rec(1,1) = pia(1)
       br1_rec(2,2) = pia(2)
       br1_rec(3,3) = pia(3)

       br2_rec(1,:) = pia(1) * (/-1, 1, 1 /)
       br2_rec(2,:) = pia(2) * (/ 1,-1, 1 /)
       br2_rec(3,:) = pia(3) * (/ 1, 1,-1 /)

       rvfac = 4
       stru%ortho = .true.
    case ('B')
       br1_rec(1,1) = pia(1)
       br1_rec(2,2) = pia(2)
       br1_rec(3,3) = pia(3)

       br2_rec(1,:) = pia(1) * (/ 0, 1, 1 /)
       br2_rec(2,:) = pia(2) * (/ 1, 0, 1 /)
       br2_rec(3,:) = pia(3) * (/ 1, 1, 0 /)

       rvfac = 2
       stru%ortho = .true.

    case ('R')
       br1_rec(1, :) = pia(1) * (/  1, 1, -2 /)/SQ3
       br1_rec(2, :) = pia(2) * (/ -1, 1,  0 /)
       br1_rec(3, :) = pia(3)

       br2_rec = br1_rec

       rvfac = 6/SQ3
       stru%ortho = .false.

    case ('C')
       ! “defaults”, to be changed in nonorthogonal XZ-case
       br1_rec(1,1) = pia(1)
       br1_rec(2,2) = pia(2)
       br1_rec(3,3) = pia(3)

       rvfac = 2
       stru%ortho = .true.

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

    call inverse3x3(br1_rec, br1_dir); br1_dir = br1_dir*2*PI
    call inverse3x3(br2_rec, br2_dir); br2_dir = br2_dir*2*PI

    stru%vol = product(stru%a) / rvfac
    stru%gbas    = br2_rec
    stru%lat2car = br1_dir
    stru%brlat   = br2_dir

    ! make_ksym() needs gbas and friends
    call read_rsym()
    call make_ksym()
  contains
    subroutine read_rsym()
!!! This is from SRC_kgen (main.f, addinv.f)
      integer :: nsym, i, i1, i2

      read(lun, '(I4)') nsym
      allocate(stru%rsym(3,3,nsym))
      do i=1,nsym
         read(lun, '(3(3I2,/), I8)') ((stru%rsym(i1,i2,i), i1=1,3), i2=1,3)
      end do
    end subroutine read_rsym

    subroutine make_ksym()
!!! This is from SRC_kgen/sdefl.f

      real(DPk), dimension(3,3) :: a
      logical                   :: addinv
      integer                   :: i, nrsym, nksym

      addinv = .true.
      nrsym  = size(stru%rsym, 3)
!!! Redefine symmetry operations from struct with
!!! unitary transformation  u(-1) * S * u

      do i=1,nrsym
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

      if (nrsym>24 .and. addinv) &
           call croak('nsym > 24 without inversion')

      if (addinv) then
         ! add inversion
         ! FIXME: check if we should actually do this (non-sp non-so)
         nksym = 2*nrsym
      else
         nksym =   nrsym
      end if
      allocate(stru%ksym(3,3,nksym))

                  stru%ksym(:,:,       1:nrsym) =  stru%rsym(:,:, 1:nrsym)
      if (addinv) stru%ksym(:,:, nrsym+1:nksym) = -stru%rsym(:,:, 1:nrsym)

      if (stru%ortho .or. stru%lattic=='CXZ') then
         ! FIXME: why the check for CXZ??  sdefl.f has it …
         do i=1,nksym
            a = matmul(matmul(stru%brlat, stru%ksym(:,:,i)), &
                 &     stru%gbas)/2/PI
            stru%ksym(:,:,i) = nint(a)
         end do
      end if
    end subroutine make_ksym
  end subroutine struct_read_unit

  subroutine struct_read_fname(fname, stru)
    use util, only: newunit

    character(*), intent(in)  :: fname
    type(struct), intent(out) :: stru

    integer :: lun

    open(unit=newunit(lun), file=fname, status='old')
    call struct_read_unit(lun, stru)
    close(lun)
  end subroutine struct_read_fname
end module structmod


!---------------- Reading ‘klist’ files            --------------------------
module kpoints
  use const,     only: BUFSZ, DPk
  use structmod, only: struct
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
contains
!!! Read standard ‘klist’
  subroutine get_kmesh_klist_unit(unit, kpoints, stru, mp_grid)
    integer,      intent(in)               :: unit
    real(DPk),    intent(out), allocatable :: kpoints(:,:)
    type(struct), intent(in)               :: stru
    integer,      intent(out), optional    :: mp_grid(3)

    integer          :: nk, ik
    character(BUFSZ) :: buf
    integer          :: kfrac(3), div
    real(DPk)        :: RR(3,3)

    if (stru%ortho) then
       RR(1,:) = stru%brlat(1,:) / stru%a
       RR(2,:) = stru%brlat(2,:) / stru%a
       RR(3,:) = stru%brlat(3,:) / stru%a
    end if

    nk = count_kmesh_klist(unit)
    allocate(kpoints(nk, 3))

    ik = 0
    readk: do
       read(unit, '(A)') buf
       if (trim(buf) == "END") exit

       ik = ik+1

       if (ik == 1 .and. present(mp_grid)) then
          read(buf, '(T86, 3I3)') mp_grid
       end if

       read(buf(11:), *) kfrac, div ! skip running number
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
    type(struct),           intent(in)  :: stru
    character, allocatable, intent(out) :: knames(:)

    integer          :: nk, ik
    character(BUFSZ) :: buf
    integer          :: kfrac(3), div
    real(DPk)        :: RR(3,3)

    if (stru%ortho) then
       RR(1,:) = stru%brlat(1,:) / stru%a
       RR(2,:) = stru%brlat(2,:) / stru%a
       RR(3,:) = stru%brlat(3,:) / stru%a
    end if

    nk = count_kmesh_band(unit)

    allocate(kpoints(nk, 3), knames(nk))

    ik  = 0
    readk: do
       read(unit, '(A)') buf
       if (trim(buf) == "END") exit
       if (buf(1:6)  == "")    cycle

       ! now we have a k-vector we should save
       ik = ik+1
       knames(ik) = buf(1:1)
       read(buf(7:), *) kfrac, div
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
       if (trim(buf) == "END") exit
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
       if (trim(buf) == "END") exit
       ! we only count k-points with a label, assuming that these
       ! delimit "paths" that will be passed to wannier90
       if (buf /= "") nk = nk+1
    end do countk
    rewind(unit)
  end function count_kmesh_band_unit

!!! Wrappers for calling the above with a filename instead of an open unit
  subroutine get_kmesh_klist_fname(fname, kpoints, stru, mp_grid)
    character(*), intent(in)               :: fname
    real(DPk),    intent(out), allocatable :: kpoints(:,:)
    type(struct), intent(in)               :: stru
    integer,      intent(out), optional    :: mp_grid(3)

    integer :: lun

    open(unit=newunit(lun), file=fname, status='old')
    call get_kmesh_klist_unit(lun, kpoints, stru, mp_grid)
    close(lun)
  end subroutine get_kmesh_klist_fname

  subroutine get_kmesh_band_fname(fname, kpoints, stru, knames)
    character(*),           intent(in)  :: fname
    real(DPk), allocatable, intent(out) :: kpoints(:,:)
    type(struct),           intent(in)  :: stru
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
end module kpoints


!---------------- Reading ‘chk’ files                  ----------------------
!
! Mostly copied from Wannier90
module Wannier90
  use const, only: DPk
  use util,  only: newunit

  implicit none
  private
  public :: chk_file, chk_read

  interface chk_read
     module procedure chk_read_fname, chk_read_unit, chk_read_argstr
  end interface chk_read

  type chk_file
     character(len=33)          :: header ! usually, date and time
     character(len=20)          :: checkpoint
     !                             #k,       #WF       #B
     integer                    :: num_kpts, num_wann, num_bands
     ! total number of neighbors per k-point
     integer                    :: nntot

     ! selection of the optimal subspace (#B × #WF × #k)
     complex(DPk), allocatable  :: u_matrix_opt(:,:,:)
     ! rotation to the optimally smooth states in subspace (#WF × #WF × #k)
     complex(DPk), allocatable  :: u_matrix    (:,:,:)
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
  end type chk_file

contains
  subroutine chk_read_unit(lun, chk, read_mmn)
    integer,           intent(in)  :: lun
    type(chk_file),    intent(out) :: chk
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
    type(chk_file),    intent(out) :: chk
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
    type(chk_file),    intent(out) :: chk
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


!---------------- Things specific to woptic            ----------------------
module woptic
  use const,     only: DPk

  implicit none
  private
  public :: get_mommat, get_mommat_1k, inwop, inwop_read

  interface get_mommat
     module procedure get_mommat_unit, get_mommat_fname
  end interface get_mommat

  interface inwop_read
     module procedure inwop_read_fname, inwop_read_unit, inwop_read_argstr
  end interface inwop_read

!!! Type representing ‘inwop’ files
  type inwop
     character(4)         :: mode
     integer              :: matelmode
     logical              :: intrahop
     real(DPk)            :: Emax, dE, deltino, beta
     integer              :: NEmin, NEmax, NEmin_w2k, NEmax_w2k
     real(DPk)            :: chempot
     real(DPk)            :: drudesep
     logical              :: orbresolv
     logical              :: selfE
     integer, allocatable :: iself(:)
     logical              :: wfrot
     real(DPk)            :: Eshift
     integer, allocatable :: ishift(:)
     ! “derived” options
     logical              :: optcond, dos, joint, peierls, need_umatrix
     logical              :: read_mommat, read_vk, read_energy, read_ham
  end type inwop

  integer, public, parameter :: &
       MODE_PEIERLS=1, MODE_INTERP=2, MODE_FULL=3, MODE_BLOCH=4, MODE_LDA=5
contains
  subroutine get_mommat_unit(unit, Vx, Vy, Vz, nk, nbmin, nbmax)
    integer,      intent(in)                    :: unit, nk, nbmin, nbmax
    complex(DPk), intent(out), dimension(:,:,:) :: Vx, Vy, Vz

    integer :: jk

    read(unit,*)

    do jk=1,nk
       call get_mommat_1k(unit, &
            &             Vx(:,:,jk), Vy(:,:,jk), Vz(:,:,jk), &
            &             nbmin, nbmax)
    end do
  end subroutine get_mommat_unit

  subroutine get_mommat_1k(unit, Vkx, Vky, Vkz, nbmin, nbmax)
    integer,      intent(in)                  :: unit, nbmin, nbmax
    complex(DPk), intent(out), dimension(:,:) :: Vkx, Vky, Vkz

    integer   :: nbminloc, nbmaxloc, i, j, nb

    nb = nbmax - nbmin + 1

    read(unit,*)
    read(unit,"(T28, 2I5)") nbminloc, nbmaxloc
    read(unit,*)

    ! skip leading junk
    do i = 1, &               ! "little Gauss":
         & (nbmin - nbminloc)*(2*nbmaxloc-nbmin-nbminloc+3)/2
       read(unit,*)
    end do

    ! here comes the juicy part
    do i=nbmin,nbmax
       do j=i,nbmaxloc
          ! skip intermediate junk
          if ((j.lt.nbmin) .or. (j.gt.nbmax)) then
             read(unit, *)
             cycle
          end if

          ! read what we need
          read(unit, '(11X,6E13.6)') Vkx(i-nbmin+1, j-nbmin+1), &
               &                     Vky(i-nbmin+1, j-nbmin+1), &
               &                     Vkz(i-nbmin+1, j-nbmin+1)
       end do
    end do

    ! skip trailing junk up to next k-point
    do i = 1, &
         & (nbmaxloc-nbmax)*(nbmaxloc-nbmax+1)/2
       read(unit,*)
    end do

    !! Complete V by hermiticity
    forall (i=1:nb, j=1:nb, j>i)
       Vkx(j, i) = conjg(Vkx(i, j))
       Vky(j, i) = conjg(Vky(i, j))
       Vkz(j, i) = conjg(Vkz(i, j))
    end forall
  end subroutine get_mommat_1k

!!! Wrapper for calling the above with a filename instead of an open
!!! unit
  subroutine get_mommat_fname(fname, Vx, Vy, Vz, nk, nbmin, nbmax)
    use util, only: newunit

    integer,      intent(in)                    :: nk, nbmin, nbmax
    complex(DPk), intent(out), dimension(:,:,:) :: Vx, Vy, Vz
    character(*), intent(in)                    :: fname

    integer :: lun

    open(unit=newunit(lun), file=fname, status='old')
    call get_mommat_unit(lun, Vx, Vy, Vz, nk, nbmin, nbmax)
    close(lun)
  end subroutine get_mommat_fname


!!! Read and check ‘inwop’ file
  subroutine inwop_read_unit(lun, inw)
    use clio, only: croak
    use util, only: uppercase, string

    integer,     intent(in)  :: lun
    type(inwop), intent(out) :: inw

    integer :: Nshift, Nself, Nb
    logical :: scissors, mixed

    read(lun,*) inw%mode, inw%matelmode, inw%intrahop
    read(lun,*) inw%Emax, inw%dE, inw%deltino, inw%beta
    read(lun,*) inw%NEmin_w2k, inw%NEmax_w2k, inw%NEmin, inw%NEmax
    read(lun,*) inw%chempot
    read(lun,*) inw%drudesep, inw%orbresolv

    read(lun,*) inw%selfE, Nself
    if (.not. inw%selfE) Nself=0
    inw%selfE = Nself/=0
    allocate(inw%iself(Nself))
    read(lun,*) inw%iself

    read(lun,*) inw%wfrot, scissors, inw%Eshift, Nshift
    if (.not. scissors .or. inw%Eshift==0) Nshift=0
    allocate(inw%ishift(Nshift))
    if (Nshift>0) read(lun,*) inw%ishift

    !! consistency checks and derived options
    inw%optcond = .true.
    inw%dos     = .true.
    inw%joint   = .false.
    select case (uppercase(inw%mode))
    case('OPT') ;
    case('DOS') ; inw%optcond = .false.
    case('JOIN'); inw%joint   = .true.
    case default; call croak("unknown MODE: "//inw%mode)
    end select

    mixed = inw%NEmin/=inw%NEmin_w2k .or. inw%NEmax/=inw%NEmax_w2k

    inw%read_energy  = mixed .or. inw%matelmode >= MODE_LDA
    inw%read_mommat  = inw%matelmode >= MODE_FULL .or. &
             &        (inw%matelmode == MODE_INTERP .and. mixed)
    inw%read_vk      = inw%matelmode == MODE_INTERP
    inw%read_ham     = inw%matelmode <  MODE_LDA
    inw%need_umatrix = inw%matelmode >  MODE_INTERP .and. &
             &         inw%matelmode <  MODE_LDA
    inw%peierls      = inw%matelmode == MODE_PEIERLS

    if (Nshift==0) inw%Eshift = 0

    if (inw%matelmode /= 1) inw%intrahop = .false.

    if (inw%matelmode < MODE_PEIERLS .or. inw%matelmode > MODE_LDA) &
         call croak("unknown matelmode: "//trim(string(inw%matelmode)))
    if (inw%Emax <= 0) &
         call croak("Emax should be >0, not " // trim(string(inw%Emax)))
    if (inw%dE <= 0) &
         call croak("dE should be >0, not " // trim(string(inw%dE)))
    if (inw%beta < 0) &
         call croak("beta should be ≥0, not " // trim(string(inw%beta)))
    if (inw%deltino <= 0) &
         call croak("deltino should be >0, not" // trim(string(inw%deltino)))
    if ( inw%Nemin_w2k<=0            .or. &
         inw%Nemin<inw%Nemin_w2k     .or. inw%NEmax<inw%NEmin .or. &
         inw%Nemax_w2k<inw%Nemin_w2k .or. inw%Nemax>inw%Nemax_w2k) &
         call croak("1 ≤ Nemin_w2k ≤ Nemin ≤ Nemax ≤ Nemax_w2k must hold")
    if (inw%drudesep < 0) &
         call croak("drudesep should be ≥0, not "//trim(string(inw%drudesep)))
    Nb = inw%Nemax_w2k - inw%Nemin_w2k + 1
    if (any(inw%ishift<= 0)) call croak("shift indices must be ≥ 1")
    if (any(inw%iself <= 0)) call croak("interacting indices must be ≥ 1")
    if (any(inw%ishift> Nb)) call croak("shift indices must be ≤ #bands")
    if (any(inw%iself > Nb)) call croak("interacting indices must be ≤ #bands")
    if (Nshift > Nb) call croak("more shift indices than bands given")
    if (Nself  > Nb) call croak("more interacting indices than bands given")
  end subroutine inwop_read_unit

!!! Wrappers for calling the above with a filename instead of an open
!!! unit
  subroutine inwop_read_fname(fname, inw)
    use util, only: newunit

    character(*), intent(in)  :: fname
    type(inwop),  intent(out) :: inw

    integer :: lun

    open(unit=newunit(lun), file=fname, status='old')
    call inwop_read_unit(lun, inw)
    close(lun)
  end subroutine inwop_read_fname

  subroutine inwop_read_argstr(arg, inw)
    use util, only: newunit
    use clio, only: argstr

    type(argstr), intent(in)  :: arg
    type(inwop),  intent(out) :: inw

    call inwop_read_fname(arg%s, inw)
  end subroutine inwop_read_argstr
end module woptic

!--------------- Files and units for woptic          ------------------------
module woptic_io
  use const, only: BUFSZ
  implicit none
  public

  integer,          parameter   :: unit_inwop=10
  character(*),     parameter   ::  suf_inwop='.inwop'
  character(BUFSZ)              ::   fn_inwop
  integer,          parameter   :: unit_fermi=11
  character(*),     parameter   ::  suf_fermi='.fermi'
  character(BUFSZ)              ::   fn_fermi
  integer,          parameter   :: unit_ham=12
  character(*),     parameter   ::  suf_ham='.hk'
  character(BUFSZ)              ::   fn_ham
  integer,          parameter   :: unit_intrahop=13
  character(*),     parameter   ::  suf_intrahop='.intrahop'
  character(BUFSZ)              ::   fn_intrahop
  integer,          parameter   :: unit_wfrot=14
  character(*),     parameter   ::  suf_wfrot='.wfrot'
  character(BUFSZ)              ::   fn_wfrot
  integer,          parameter   :: unit_hr=15
  character(*),     parameter   ::  suf_hr='_hr.dat'
  character(BUFSZ)              ::   fn_hr
  integer,          parameter   :: unit_chk=16
  character(*),     parameter   ::  suf_chk='.chk'
  character(BUFSZ)              ::   fn_chk
  integer,          parameter   :: unit_struct=20
  character(*),     parameter   ::  suf_struct='.struct'
  character(BUFSZ)              ::   fn_struct
  integer,          parameter   :: unit_symop=21
  character(*),     parameter   ::  suf_symop='.symop'
  character(BUFSZ)              ::   fn_symop
  integer,          parameter   :: unit_energy=22
  character(*),     parameter   ::  suf_energy='.energy'
  character(BUFSZ)              ::   fn_energy
  integer,          parameter   :: unit_mommat=23
  character(*),     parameter   ::  suf_mommat='.mommat2'
  character(BUFSZ)              ::   fn_mommat
  integer,          parameter   :: unit_vk(3)=(/ 24, 25, 26 /)
  character(*),     parameter   ::  suf_vk(3)= (/ '.vkx', '.vky', '.vkz' /)
  character(BUFSZ)              ::   fn_vk(3)
  integer,          parameter   :: unit_vr=27
  character(*),     parameter   ::  suf_vr='.vr'
  character(BUFSZ)              ::   fn_vr
  integer,          parameter   :: unit_klist=30
  character(*),     parameter   ::  suf_klist='.klist'
  character(BUFSZ)              ::   fn_klist
  integer,          parameter   :: unit_fklist=31
  character(*),     parameter   ::  suf_fklist='.klist_full'
  character(BUFSZ)              ::   fn_fklist
  integer,          parameter   :: unit_kadd=32
  character(*),     parameter   ::  suf_kadd='.klist_add'
  character(BUFSZ)              ::   fn_kadd
  integer,          parameter   :: unit_tet=40
  character(*),     parameter   ::  suf_tet='.tetra'
  character(BUFSZ)              ::   fn_tet
  integer,          parameter   :: unit_ftet=41
  character(*),     parameter   ::  suf_ftet='.tetra_full'
  character(BUFSZ)              ::   fn_ftet
  integer,          parameter   :: unit_map=42
  character(*),     parameter   ::  suf_map='.map'
  character(BUFSZ)              ::   fn_map
  integer,          parameter   :: unit_voe=43
  character(*),     parameter   ::  suf_voe='.voe'
  character(BUFSZ)              ::   fn_voe
  integer,          parameter   :: unit_optcond=51
  character(*),     parameter   ::  suf_optcond='.optcondw'
  character(BUFSZ)              ::   fn_optcond
  integer,          parameter   :: unit_optorb(6)=(/511,512,513,514,515,516/)
  character(*),     parameter   ::  suf_optorb(6)= &
       (/ suf_optcond//'_orbxx', suf_optcond//'_orbxy', suf_optcond//'_orbxz',&
       &  suf_optcond//'_orbyy', suf_optcond//'_orbyz', suf_optcond//'_orbzz'/)
  character(BUFSZ)              ::   fn_optorb(6)
  integer,          parameter   :: unit_contr=52
  character(*),     parameter   ::  suf_contr='.kcontribw'
  character(BUFSZ)              ::   fn_contr
  integer,          parameter   :: unit_wdos=53
  character(*),     parameter   ::  suf_wdos='.wdos'
  character(BUFSZ)              ::   fn_wdos
  integer,          parameter   :: unit_doscontr=54
  character(*),     parameter   ::  suf_doscontr='.wdoskcontribw'
  character(BUFSZ)              ::   fn_doscontr
  integer,          parameter   :: unit_K1=55
  character(*),     parameter   ::  suf_K1='.K1w'
  character(BUFSZ)              ::   fn_K1
  integer,          parameter   :: unit_selfE=56
  character(*),     parameter   ::  suf_selfE='.selfE'
  character(BUFSZ)              ::   fn_selfE
  integer,          parameter   :: unit_outwop=60
  character(*),     parameter   ::  suf_outwop='.outputwop'
  character(BUFSZ)              ::   fn_outwop
  integer,          parameter   :: unit_outvr=61
  character(*),     parameter   ::  suf_outvr='.vrout'
  character(BUFSZ)              ::   fn_outvr
  integer,          parameter   :: unit_outref=62
  character(*),     parameter   ::  suf_outref='.outputref'
  character(BUFSZ)              ::   fn_outref
contains
  subroutine print1or3(lun, label, unit, tensor)
    use const, only: DPk
    implicit none

    integer,      intent(in) :: lun
    character(*), intent(in) :: label, unit
    real(DPk),    intent(in) :: tensor(3,3)

    character(*), parameter :: &
         fmt1 = '(A12, " [", A, "]", F12.3)', &
         fmt3 = '(A11, "(xx,yy,zz) [", A, "]", 3F12.3)'

    real(DPk) :: xx, yy, zz

    xx=tensor(1,1); yy=tensor(2,2); zz=tensor(3,3)

    if (max(abs((xx-yy)/xx), abs((xx-zz)/xx), abs((zz-yy)/xx)) > 1e-3) then
       write(lun, fmt3) label, unit, xx, yy, zz
    else
       write(lun, fmt1) label, unit, xx
    end if
  end subroutine print1or3

  subroutine set_casename(file, band)
    use clio,  only: argstr

    implicit none

    type(argstr), intent(in)           :: file
    logical,      intent(in), optional :: band

    character(5) :: suf
    integer      :: i
    suf=''
    if (present(band)) then
       if (band) suf='_band'
    end if

    fn_inwop    =trim(file%s)//trim(suf_inwop    )
    fn_fermi    =trim(file%s)//trim(suf_fermi    )
    fn_ham      =trim(file%s)//trim(suf_ham      )
    fn_intrahop =trim(file%s)//trim(suf_intrahop )
    fn_wfrot    =trim(file%s)//trim(suf_wfrot    )
    fn_hr       =trim(file%s)//trim(suf_hr       )
    fn_chk      =trim(file%s)//trim(suf_chk      )
    fn_struct   =trim(file%s)//trim(suf_struct   )
    fn_symop    =trim(file%s)//trim(suf_symop    )
    fn_energy   =trim(file%s)//trim(suf_energy   )
    fn_mommat   =trim(file%s)//trim(suf_mommat   )
    fn_klist    =trim(file%s)//trim(suf_klist    )//suf
    fn_vr       =trim(file%s)//trim(suf_vr       )
    fn_fklist   =trim(file%s)//trim(suf_fklist   )
    fn_kadd     =trim(file%s)//trim(suf_kadd     )
    fn_tet      =trim(file%s)//trim(suf_tet      )
    fn_ftet     =trim(file%s)//trim(suf_ftet     )
    fn_map      =trim(file%s)//trim(suf_map      )
    fn_voe      =trim(file%s)//trim(suf_voe      )
    fn_optcond  =trim(file%s)//trim(suf_optcond  )
    fn_contr    =trim(file%s)//trim(suf_contr    )//suf
    fn_wdos     =trim(file%s)//trim(suf_wdos     )//suf
    fn_doscontr =trim(file%s)//trim(suf_doscontr )
    fn_K1       =trim(file%s)//trim(suf_K1       )//suf
    fn_selfE    =trim(file%s)//trim(suf_selfE    )//suf
    fn_outwop   =trim(file%s)//trim(suf_outwop   )
    fn_outvr    =trim(file%s)//trim(suf_outvr    )
    fn_outref   =trim(file%s)//trim(suf_outref   )
    do i=1,size(fn_optorb)
       fn_optorb(i)=trim(file%s)//trim(suf_optorb(i))
    end do
    do i=1,size(fn_vk)
       fn_vk(i)    =trim(file%s)//trim(suf_vk(i)    )
    end do
  end subroutine set_casename
end module woptic_io


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
