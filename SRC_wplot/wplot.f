!!! wien2wannier/SRC_wplot/wplot.f
!!!
!!!    Main program ‘wplot’
!!!
!!! Copyright 2010-2012 Jan Kuneš, Philipp Wissgott
!!!           2014-2016 Elias Assmann

!!/--- Files expected in ‘def’ ---
!!  5 inwplot	  'old'		'formatted'
!!  6 outputwplot 'unknown'	'formatted'
!!  7 grid	  'unknown'	'formatted'	! mode ‘ANY’
!!  8 struct	  'old'		'formatted'
!! 10 vector	  'old'		'unformatted'
!! 18 vsp	  'old'		'formatted'
!! 21 psink	  'unknown'	'formatted'
!! 22 psiarg	  'unknown'	'formatted'
!! 31 inwf	  'unknown'	'formatted'
!! 32 chk	  'old'		'unformatted'
!!\---

program wf_r
  use const, only: BUFSZ
  use clio,  only: fetcharg, argstr
  use param, only: unit_def, unit_vector, unit_out
  use wplot, only: vecfn, idx_wann, iproc,             &
       &           unit_psink, unit_psiarg, unit_grid, &
       &           psinkfn, psiargfn, outfn, gridfn

  !! procedure includes
  use main_m

  implicit none

  type(argstr)     :: defname
  character(BUFSZ) :: fname
  character(11)    :: status,form
  integer          :: iunit

                                    call fetcharg(1, defname)
  if (command_argument_count() > 1) call fetcharg(2, idx_wann)
  if (command_argument_count() > 2) call fetcharg(3, iproc)

  open(unit_def, FILE=defname%s, STATUS='OLD')
  def: do
     read(unit_def,*,END=8001) iunit, fname, status, form

     select case (iunit)
     case (unit_vector)
        vecfn    = fname
     case (unit_psink)
        psinkfn  = fname
     case (unit_psiarg)
        psiargfn = fname
     case(unit_out)
        outfn    = fname
     case(unit_grid)
        gridfn   = fname
        open(iunit, FILE=fname, STATUS=status, FORM=form)
     case default
        open(iunit, FILE=fname, STATUS=status, FORM=form)
     end select
  end do def

8001 close(unit_def)

  !     << start wave function evaluation >>
  call main()
end program wf_r


!!/---
!! Local Variables:
!! mode: f90
!! End:
!!\---
!!
