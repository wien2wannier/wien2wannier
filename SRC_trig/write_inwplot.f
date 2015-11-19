!!! wien2wannier/SRC_trig/write_wplotin.f90
!!!
!!!    Prepares input seed.wplotin for wplot
!!!
!!! Copyright 2009-2012 Philipp Wissgott
!!!                2014 Elias Assmann
!!!
!!! $Id: write_inwplot.f 289 2014-10-10 12:18:40Z assmann $

PROGRAM write_wplotin
  use const, only: BUFSZ
  use clio,  only: argstr, fetcharg, croak

  implicit none
  integer          :: divisor, vec(3), iarg
  character(BUFSZ) :: dummy, woutfileend
  type(argstr)     :: file, arg

  character(*), parameter :: wplotinfileend=".inwplot"
  character(*), parameter :: &
       fmt_idiv = "(4(I0,1X), T45, A)", &
       fmt_A    = "(A, T45, A)", &
       fmt_grid = "(6(I0,1X), T45, A)"

  woutfileend  = ".wout"

  !command line argument read-in
  do iarg=1,command_argument_count()
     call fetcharg(iarg, arg)

     if (arg%s(1:1).eq.'-') then
        if (arg%s(2:3) == 'up' .or. arg%s(2:3) == 'dn') then     
           !for spin-polarized calc. the argendings have additional up/dn
           woutfileend = ".wout"//arg%s(2:3)
        elseif (arg%s(2:2).eq.'h') then
           write(*,*) "prepares input file `case.wplotin' for wplot"
           write(*,*) 'Usage: write_wplotin [-up/-dn] CASE'
           stop 
        else
           call croak("unknown option: " // trim(arg%s))
        endif
     else
        call fetcharg(iarg, file)
     endif
  enddo

  write(*,*)"From wannier90 output "//trim(file%s)//trim(woutfileend) &
       //"(in Angstrom):"

  !read wannier90 output
  open(UNIT=1, FILE=trim(file%s)//woutfileend, STATUS='old')

  do
     read(1, "(A)", END=101) dummy
     if (dummy(2:12) == "Final State") exit
  end do
  do
     read(1, "(A)", END=101) dummy
     print '(A)', trim(dummy)

     if (dummy(2:12) == "") exit
  end do
101 close(1)


  !write the input file
  print '(A)', "Enter origin of spatial mesh in fractions of the &
       & conv. unit vectors [n1 n2 n3 ndiv]"

  read (*,*) vec(:), divisor
  open(UNIT=1,FILE=trim(file%s)//wplotinfileend, STATUS='unknown')
  write(1, fmt_A) "3D ORTHO", "# mode O(RTHOGONAL)|N(ON-ORTHOGONAL)"
  write(1, fmt_idiv) vec(:), divisor, "# x,y,z,div of origin"

  print '(A)', "Enter endpoint on x-axis [n1 n2 n3 ndiv]"
  read (*,*) vec(:), divisor
  write(1, fmt_idiv) vec(:), divisor, "#              x-end"

  print '(A)', "Enter endpoint on y-axis [n1 n2 n3 ndiv]"
  read (*,*) vec(:), divisor
  write(1, fmt_idiv) vec(:), divisor, "#              y-end"

  print '(A)', "Enter endpoint on z-axis [n1 n2 n3 ndiv]"
  read (*,*) vec(:), divisor
  write(1, fmt_idiv) vec(:), divisor, "#              z-end"

  print '(A)',"Enter number of mesh points, [nx ny nz]"
  read (*,*) vec(:)
  write(1, fmt_grid) vec(:), 0,0,0, "# grid points and echo increments"

  write(1, fmt_A)"NO", "# DEP(HASING)|NO (POST-PROCESSING)"
  write(1, fmt_A)"WAN ANG LARGE", "# switch ANG|ATU|AU LARGE|SMALL"
  write(1, fmt_A)"1 1", "# k-point, Wannier index"
END PROGRAM write_wplotin

!!/---
!! Local Variables:
!! mode: f90
!! End:
!!\---
