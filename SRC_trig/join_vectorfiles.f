!!! wien2wannier/SRC_trig/join_vectorfiles.f
!!! 
!!!    Joins multiple WIEN2K vector files to one for further processing
!!!
!!!    Usage: join_vectorfiles [-up/-dn] [-c] <case> <numberofparallelfiles>
!!!
!!! Copyright 2010-2012 Philipp Wissgott
!!!           2013-2014 Elias Assmann
!!!
!!! $Id: join_vectorfiles.f 224 2014-05-30 14:56:29Z assmann $

PROGRAM join_vectorfiles
  use const,     only: R8, C16, BUFSZ
  use structmod, only: struct, struct_read
  use util,      only: paropen
  use clio,      only: argstr, fetcharg, croak
  use kpoints,   only: count_kmesh_klist

  implicit none

  integer, parameter :: unit_def=1
  integer, parameter :: unit_struct=20, unit_vector=10, unit_vector2=11
  integer, parameter :: unit_klist = 4, unit_energy=50, unit_energy2=51
  integer, parameter :: lmax=13, lomax=3, nloat=3

  integer          :: iproc=0
  character(BUFSZ) :: vecfn, enefn
  logical          :: complex=.false.

  real(R8),     allocatable :: EIGVAL(:), E(:,:), ELO(:,:,:)
  integer,      allocatable :: KZZ(:,:)
  real(R8),     allocatable :: Z (:,:)
  complex(C16), allocatable :: ZC(:,:)

  integer       :: nkpoints, NE, NV
  integer       :: jatom,i,j,jl,jj,jk
  real(R8)      :: SX, SY, SZ, weight, eorb_ind
  CHARACTER(3)  :: IPGR
  CHARACTER(10) :: KNAME
  type(struct)  :: stru

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

 allocate( E(LMAX,stru%nneq) )
 allocate( ELO(0:LOMAX,nloat,stru%nneq) )
 
 files: do j=1,iproc
    call paropen(unit_vector, vecfn, iproc, j, FORM='unformatted')
    call paropen(unit_energy, enefn, iproc, j, FORM=  'formatted')
 
   atoms: do jatom= 1,stru%nneq
       read(unit_vector) (E(jl,jatom),jl=1,LMAX)
       read(unit_vector) ((ELO(jl,jk,jatom),jl=0,LOMAX),jk=1,nloat)
       read(unit_energy,'(100(f9.5))') (E(jl,jatom),jl=1,LMAX),eorb_ind
       read(unit_energy,'(100(f9.5))') ((ELO(jl,jk,jatom),jl=0,LOMAX),jk=1,nloat)
       if (j.eq.1) then
          write(unit_vector2) (E(jl,jatom),jl=1,LMAX)
          write(unit_vector2) ((ELO(jl,jk,jatom),jl=0,LOMAX),jk=1,nloat)
          write(unit_energy2,'(100(f9.5))') (E(jl,jatom),jl=1,LMAX),eorb_ind
          write(unit_energy2,'(100(f9.5))') ((ELO(jl,jk,jatom),jl=0,LOMAX),jk=1,nloat)
       endif
    enddo atoms

    k_points: do jk=1,nkpoints ! /nfiles
       read(unit_vector,end=888,err=888) SX, SY, SZ, KNAME, NV, NE, WEIGHT!, IPGR
       allocate( KZZ(3,NV) )
       read(unit_energy,'(3e19.12,a10,2i6,f5.1,a3)') SX, SY, SZ, KNAME, NV, NE, WEIGHT, IPGR 
       !        write(*,*)SX, SY, SZ, KNAME, NV, NE, WEIGHT, IPGR 
       read(unit_vector) (KZZ(1,I),KZZ(2,I),KZZ(3,I),I=1,NV)

       write(unit_vector2)SX, SY, SZ, KNAME, NV, NE, WEIGHT!, IPGR
       write(unit_energy2,'(3e19.12,a10,2i6,f5.1,a3)') &
            SX, SY, SZ, KNAME, NV, NE, WEIGHT, IPGR 
       write(unit_vector2) (KZZ(1,I),KZZ(2,I),KZZ(3,I),I=1,NV)
       if (complex) then
          allocate(ZC(NV,NE))
       else
          allocate(Z (NV,NE))
       end if
       allocate(EIGVAL(NE))

       do jj = 1, NE
          read(unit_vector) I, EIGVAL(I)
          write(unit_vector2)I,EIGVAL(I)

          read(unit_energy,*) I, EIGVAL(I)
          write(unit_energy2,*)I,EIGVAL(I)

          if (complex) then
             read(unit_vector) (ZC(jl,I),jl=1,NV)
             write(unit_vector2) (ZC(jl,I),jl=1,NV)
          else
             read(unit_vector) (Z(jl,I),jl=1,NV)
             write(unit_vector2) (Z(jl,I),jl=1,NV)
          endif
       enddo

       deallocate(EIGVAL,KZZ)
       if (complex) then
          deallocate(ZC)
       else
          deallocate(Z)
       end if
888    continue
    enddo k_points
    close(unit_vector)
    close(unit_energy)
 enddo files
END PROGRAM

!!/---
!! Local Variables:
!! mode: f90
!! End:
!!\---
