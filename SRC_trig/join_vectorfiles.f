!!! wien2wannier/SRC_trig/join_vectorfiles.f
!!!
!!!    Joins multiple WIEN2K vector files to one for further processing
!!!
!!!    Usage: join_vectorfiles [-up/-dn] [-c] <case> <numberofparallelfiles>
!!!
!!! Copyright 2010-2012 Philipp Wissgott
!!!           2013-2016 Elias Assmann

program join_vectorfiles
  use const,     only: R8, C16, BUFSZ
  use structmod, only: struct_t, struct_read
  use util,      only: paropen
  use clio,      only: argstr, fetcharg, croak
  use kpoints,   only: count_kmesh_klist

  implicit none

  integer, parameter :: unit_def=1
  integer, parameter :: unit_struct=20, unit_vector=10, unit_vector2=11
  integer, parameter :: unit_klist = 4, unit_energy=50, unit_energy2=51
  integer, parameter :: lmax=13, lomax=3, nloat=3

  character(*), parameter :: fmt_header = '(100(f9.5))'
  character(*), parameter :: fmt_kpt    = '(3e19.12,a10,2i6,f5.1,a3)'

  integer          :: iproc=0
  character(BUFSZ) :: vecfn='', enefn=''
  logical          :: complex=.false., do_vector=.false.

  real(R8),     allocatable :: E(:), ELO(:,:)
  integer,      allocatable :: KZZ(:,:)
  real(R8),     allocatable :: Z  (:)
  complex(C16), allocatable :: ZC (:)

  integer        :: nkpoints, NE, NV
  integer        :: jatom,i,j,jj,jk
  real(R8)       :: SX, SY, SZ, weight, eorb_ind, eigval
  CHARACTER(3)   :: IPGR
  CHARACTER(10)  :: KNAME
  type(struct_t) :: stru

  type(argstr)     :: defname, cmplxarg
  character(BUFSZ) :: fname
  character(11)    :: status, form
  integer          :: iunit

  call fetcharg(1, defname)
  if (command_argument_count() > 1) then
     call fetcharg(2, iproc)
  end if
  if (command_argument_count() > 2) then
     call fetcharg(3, cmplxarg)
     if (cmplxarg%s == '-c') then
        complex = .true.
     else
        call croak("join_vectorfiles: illegal argument: " // trim(cmplxarg%s))
     end if
  end if

  open(unit_def, FILE=defname%s, STATUS='OLD')
  def: do
     read(unit_def,*,END=8001) iunit, fname, status, form

     select case (iunit)
     case (unit_vector)
        vecfn = fname
        do_vector = .true.
     case (unit_energy)
        enefn = fname
     case default
        open(iunit, FILE=fname, STATUS=status, FORM=form)
     end select
  end do def

8001 close(unit_def)

 call struct_read(unit_struct, stru)
 close(unit_struct)

 nkpoints = count_kmesh_klist(unit_klist)
 close(unit_klist)

 allocate( E  (LMAX) )
 allocate( ELO(0:LOMAX,nloat) )

 files_ene: do j=1,iproc
    call paropen(unit_energy, enefn, iproc, j, FORM='formatted')

   atoms_ene: do jatom= 1,stru%nneq
       read(unit_energy, fmt_header) E, eorb_ind
       read(unit_energy, fmt_header) ELO
       if (j.eq.1) then
          write(unit_energy2, fmt_header) E, eorb_ind
          write(unit_energy2, fmt_header) ELO
       endif
    enddo atoms_ene

    k_points_ene: do jk=1,nkpoints
       read (unit_energy,  fmt_kpt, end=101) &
            SX, SY, SZ, KNAME, NV, NE, WEIGHT, IPGR
       write(unit_energy2, fmt_kpt) &
            SX, SY, SZ, KNAME, NV, NE, WEIGHT, IPGR

       eig_ene: do jj = 1, NE
          read (unit_energy, *) I, EIGVAL
          write(unit_energy2,*) I, EIGVAL
       end do eig_ene
    end do k_points_ene
101 close(unit_energy)
 end do files_ene

 do_vec: if (do_vector) then
    files_vec: do j=1,iproc
       call paropen(unit_vector, vecfn, iproc, j, FORM='unformatted')

       atoms_vec: do jatom= 1,stru%nneq
          read(unit_vector) E
          read(unit_vector) ELO

          if (j.eq.1) then
             write(unit_vector2) E
             write(unit_vector2) ELO
          endif
       end do atoms_vec

       k_points_vec: do jk=1,nkpoints
          read (unit_vector, end=102)  SX, SY, SZ, KNAME, NV, NE, WEIGHT
          write(unit_vector2)          SX, SY, SZ, KNAME, NV, NE, WEIGHT

          allocate( KZZ(3,NV) )

          if (complex) then
             allocate(ZC(NV))
          else
             allocate(Z (NV))
          end if

          read (unit_vector)  KZZ
          write(unit_vector2) KZZ

          do jj = 1, NE
             read (unit_vector)  I, EIGVAL
             write(unit_vector2) I, EIGVAL

             if (complex) then
                read  (unit_vector)  ZC
                write (unit_vector2) ZC
             else
                read (unit_vector)  Z
                write(unit_vector2) Z
             endif
          end do

          deallocate(KZZ)
          if (complex) then
             deallocate(ZC)
          else
             deallocate(Z)
          end if
       end do k_points_vec
102    close(unit_vector)
    end do files_vec
 end if do_vec
end program join_vectorfiles


!!/---
!! Local Variables:
!! mode: f90
!! End:
!!\---
