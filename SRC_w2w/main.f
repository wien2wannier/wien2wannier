!!! wien2wannier/SRC_w2w/main.f
!!!
!!!    Main program ‘w2w’
!!!
!!! Copyright 2010-2012 Jan Kuneš, Philipp Wissgott
!!!           2013-2016 Elias Assmann

!!/--- Files expected in ‘def’ ---
!!  5 inwf	'old'	  'formatted'
!!  6 outputwf	'unknown' 'formatted'
!!  7 amn	'unknown' 'formatted'
!!  8 mmn	'unknown' 'formatted'
!!  9 vector	'unknown' 'unformatted'
!! 10 nnkp	'old'	  'formatted'
!! 12 eig	'unknown' 'formatted'
!! 18 vsp	'old'	  'formatted'
!! 20 struct	'old'	  'formatted'
!! 50 energy	'old'     'formatted'
!! 51 fermi	'old'     'formatted'
!!\---

program wf
  use param,     only: unit_def, unit_vector, unit_in, unit_out, unit_struct, &
       &               wien2wannier_version, Lmax2, Nrad, Nrf
  use w2w,       only: Nmat, unit_nnkp, unit_amn, unit_mmn, unit_ene, &
       &               unit_fermi
  use const,     only: R8, BUFSZ
  use xa,        only: init_xa
  use xa3,       only: init_xa3
  use bessel,    only: init_bessel
  use Amn_Mmn,   only: init_Amn_Mmn
  use pairs,     only: kp, kpb, bqx,bqy,bqz, bqx1,bqy1,bqz1, init_pairs
  use util,      only: paropen, ptime
  use wien2k,    only: errflg, errclr, gtfnam
  use clio,      only: croak
  use structmod, only: struct_t, struct_read
  use inwfmod,   only: inwf_t, inwf_read

  !! procedure includes
  use read_vec_m
  use gaunt2_m
  use planew_m
  use l2mmn_m
  use l2amn_m

  implicit none

  type(inwf_t) :: inwf

  character(len=   11)  :: status,form
  character(len=BUFSZ)  :: deffn, errfn, aline
  character(len=BUFSZ)  :: fname, vecfn, enefn, iomsg

  integer :: iloop, i, j, ios, iproc, irecl, iunit
  integer :: maxx,maxy,maxz, maxg, n, n_pair, nen
  integer :: Nb, Nk, nntot, kkk, iostat

  real(r8) :: efermi

  type(struct_t) :: stru

!-----------------------------------------------------------------------

  call gtfnam(deffn,errfn,iproc)
  call errflg(errfn,'Error in W2W')
  open(unit_def, FILE=deffn, STATUS='old')
  def: do
     read(unit_def, *, end=20) iunit, fname, status, form, irecl

     select case (iunit)
     case (unit_vector)
        vecfn=fname
     case (50)
        enefn=fname
     case default
        open(iunit, FILE=fname, STATUS=status, FORM=form, &
             IOSTAT=iostat, IOMSG=iomsg)

        if(iostat /= 0) call croak('error while processing def file `' &
             // trim(deffn) // "': " // trim(iomsg))
     end select
  end do def
20 close(unit_def)

  write(unit_out, '("W2W ", A /)') wien2wannier_version

!!!.....READ STRUCT
  call struct_read(unit_struct, stru)

!!!....Find nmat and Nk
!!! Nk could be gotten easier from ‘klist’ -- can we get nmat
!!! somewhere else?  (‘vector’?)
  Nk=0

  enefile: do iloop=1,max(iproc,1)
     call paropen(unit_ene, ENEFN, iproc, iloop, STATUS='old')

     header: do I=1,stru%Nneq
        read(unit_ene,*)
        read(unit_ene,*)
     end do header
     kpts: do
        read(unit_ene,'(67X, 2i6)', IOSTAT=ios) N, NEn
        if (ios /= 0) exit kpts
        Nk=Nk+1
        nmat=max(n,nmat)
        do i=1,NEn
           read(unit_ene,*)
        enddo
     end do kpts

     close(unit_ene)
  end do enefile

!!!.....READ INPUT AND POTE
  nnkpt: do
     read(unit_nnkp,'(a80)',end=114) aline
     if(index(aline,'begin kpoints').ne.0) then
        read(unit_nnkp,*) n
        if (n/=Nk) call croak('inconsistent numbers of k-points between&
             & nnkp and energy files')
        cycle nnkpt
     endif

     if(index(aline,'begin nnkpts').ne.0) then
        read(unit_nnkp,*) nntot
        n_pair=Nk*nntot
        call init_pairs(n_pair)
        do i=1,N_pair
           read(unit_nnkp,*) kp(i),kpb(i),bqx1(i),bqy1(i),bqz1(i)
        enddo
        exit nnkpt
     endif
  end do nnkpt
114 close(unit_nnkp)

  write(unit_out,*)'NUM_KPTS=',Nk
  write(unit_out,*)'NNTOT=',NNTOT
  write(unit_out,*)'N_pair=',N_pair

  call inwf_read(unit_in, inwf)
  Nb = inwf%bmax - inwf%bmin + 1

  write(unit_out, "(' MODE=')", ADVANCE='no')
  if (inwf%Mmn) write(unit_out, '(" Mmn")', ADVANCE='no')
  if (inwf%Amn) write(unit_out, '(" Amn")', ADVANCE='no')
  write(unit_out,*)
  write(unit_out, '(" band window = [", I0, ", ", I0, "]")') &
       inwf%bmin, inwf%bmax

  call init_xa3(Nb,nmat,Nk+1)
  call init_xa(LMAX2,NMAT,NRAD,Nb)
  call init_Amn_Mmn(Nb, N_pair)

  read(unit_fermi, *) efermi

  kkk=0
  maxx=0; maxy=0; maxz=0
  vectorfiles: do iloop=1, max(iproc, 1)
     call paropen(unit_vector, vecfn, iproc, iloop, &
          STATUS='old', FORM='unformatted')

     call read_vec(inwf%Bmin, inwf%Bmax, stru%Nneq, &
          &        kkk, maxx, maxy, maxz, Efermi)

     close(unit_vector)
  end do vectorfiles

  if(kkk /= Nk) call croak('inconsistent numbers of k-points between&
       & vector and energy files')

  read_proj: if (inwf%Amn) then
     write(unit_out, *)
     write(unit_out,*)'Initial orbital projections:'

     do i=1,inwf%Nproj
        write(unit_out, *)
        if (i==1) then
           write(unit_out, "(' atom  l  m    Re      Im')")
           write(unit_out, "(' ', 28('-'))")
        end if

        do j=1, inwf%projections(i)%NY
           write(unit_out, '(i5,2x,i1,1x,i2,"  (",f6.3,", ",f6.3,")"  )') &
                inwf%projections(i)%iat(j), inwf%projections(i)%l    (j), &
                inwf%projections(i)%m  (j), inwf%projections(i)%coeff(j)
        enddo
     enddo
  endif read_proj

  call init_bessel(Lmax2, inwf%LJmax, Nrad, Nrf)
  call gaunt2

  write(unit_out, "(////,30X,50(1H-),/,33X,'S T R U C T U R A L   ', &
         & 'I N F O R M A T I O N',/,30X,50(1H-),//)")
  write(unit_out, "(3X,'SUBSTANCE',20X,'= ',A80,/)") stru%title

  if (inwf%Amn) then
     write(unit_amn,'(A20)')  stru%title
     write(unit_amn,'(3I12)') Nb, Nk, inwf%Nproj
  endif
  if (inwf%Mmn) then
     write(unit_mmn,'(A20)')  stru%title
     write(unit_mmn,'(3I12)') Nb, Nk, Nntot
  endif

  write(unit_out, "(3X,'LATTICE',22X,'= ',A4)")                  stru%lattic
  write(unit_out, "(3X,'LATTICE CONSTANTS ARE',8X,'= ',3F12.7)") stru%a
  write(unit_out, "(3X,'NUMBER OF ATOMS IN UNITCELL  = ',I3)")   stru%Nneq
  write(unit_out, "(3X,'MODE OF CALCULATION IS',7X,'= ',A4)")    stru%mode

  write(unit_out,*)' gbas,  stru2frac'
  do i=1,3
     write(unit_out, '(3f10.5,3x,3i3)') &
          stru%conv_rec(i,:), stru%stru2frac(i,:)
  end do

  !     rotate boundary vectors
  do i=1,N_pair
     bqx(i) = dot_product(stru%stru2frac(1,:), (/bqx1(i), bqy1(i), bqz1(i)/))
     bqy(i) = dot_product(stru%stru2frac(2,:), (/bqx1(i), bqy1(i), bqz1(i)/))
     bqz(i) = dot_product(stru%stru2frac(3,:), (/bqx1(i), bqy1(i), bqz1(i)/))
  enddo

  !.....CALCULATE CHARGE DENSITY CLM(R) IN SPHERES,  PARTIAL CHARGES

  ! l2mmn, l2amn need vector file to be open
  call paropen(unit_vector, vecfn, iproc, 1, &
       &       STATUS='old', FORM='unformatted')

  call ptime(unit_out)

  if (inwf%Mmn) then
     call l2mmn(stru, inwf, Nk, Nntot)
     call ptime('l2Mmn')
     !JXZ: MAXG is not pre-defined
     MAXG = 0
     write(unit_out,*)'MXG=',MAXG
     call planew(stru, Nb, Nk, Nntot, maxx+1, maxy+1, maxz+1)
     call ptime('planew')
  endif

  if (inwf%Amn) then
     call ptime(unit_out)
     call l2amn(stru, inwf, Nk)
     call ptime('l2Amn')
  endif

  call ERRCLR(ERRFN)
  print "('W2W END')"
end program wf


!!/---
!! Local Variables:
!! mode: f90
!! End:
!!\---
!!
