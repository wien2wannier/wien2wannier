!!! wien2wannier/SRC_trig/convham.f
!!!
!!!    Program to Fourier transform the real space Hamiltonian given
!!!    in case_hr.dat to k space H(k) where the k-list is given by
!!!    case.klist
!!!
!!! Copyright 2009-2012 Philipp Wissgott
!!!           2014-2015 Elias Assmann

program convert_hamiltonian
  use wien2k,    only: gtfnam, errflg, errclr
  use kpoints,   only: get_kmesh_klist
  use clio,      only: croak
  use const,     only: DPk, TAU, BUFSZ
  use structmod, only: struct_t, struct_read

  implicit none

  integer, parameter :: unit_klist=4, unit_struct=8, unit_hr=9, unit_hk=10

  character(len=*), parameter ::                         &
       fmt_hk_head = '(I10,2I6," #kpoints wann bands")', &
       fmt_hk_num  = 'E21.12'

  type(struct_t) :: stru

  real(DPk),    allocatable :: kpts(:,:), rweights(:), rvec(:,:)
  complex(DPk), allocatable :: Hr(:,:,:), Hk(:,:)
  complex(DPk)              :: ee

  character(BUFSZ) :: deffn, errfn, status, form, fname, hkfmt

  real(DPk) :: x, y, rdotk
  integer   :: num_wann, nrpts
  integer   :: iunit, ir, iw, jw, ik, i, j

  call gtfnam(deffn, errfn, i)
  call errflg(errfn, 'Error in CONVERT_HAMILTONIAN')
  open (1, FILE=deffn, STATUS='old')
  do
     read(1, *, END=110) iunit, fname, status, form
     open(iunit, FILE=fname, STATUS=status, FORM=form)
  end do
110 continue

  call struct_read(unit_struct, stru)
  call get_kmesh_klist(unit_klist, kpts, stru)

  read(unit_hr,*)
  read(unit_hr,*) num_wann
  read(unit_hr,*) nrpts
  allocate(rweights(nrpts), rvec(nrpts,3),&
       &   Hr(nrpts, num_wann, num_wann), Hk(num_wann, num_wann))
  read(unit_hr,*) rweights
  do ir=1,nrpts
     do jw=1,num_wann
        do iw=1,num_wann
           read(unit_hr,*) rvec(ir,:), i, j, x, y
           Hr(ir, iw, jw) = cmplx(x, y, DPk)

           if (i/=iw .or. j/=jw) &
                call croak('inconsistency in reading _hr.dat')
        end do
     end do
  end do

  write(unit_hk, fmt_hk_head) &
       size(kpts, 1), num_wann, num_wann

  write(hkfmt, '("(", I0, A, ")")') 2*num_wann, fmt_hk_num

  do ik = 1,size(kpts,1)
     Hk = 0
     do ir = 1,nrpts
        rdotk = TAU*dot_product(kpts(ik,:), rvec(ir,:))
        ee = exp((0,1)*rdotk) / rweights(ir)

        Hk(:,:) = Hk(:,:) + ee*Hr(ir,:,:)
     end do

     write(unit_hk, '(3F12.8)') kpts(ik,:)
     do iw=1,num_wann
        write(unit_hk, hkfmt) Hk(iw,:)
     end do
  end do

  call errclr(errfn)
end program convert_hamiltonian


!!/---
!! Local Variables:
!! mode: f90
!! End:
!!\---
!!
!! Time-stamp: <2016-07-06 16:18:50 assman@faepop71.tu-graz.ac.at>
