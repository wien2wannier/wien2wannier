!!! wien2wannier/SRC_trig/shifteig.f
!!!
!!!    Shifts the energy within the seed.eig file to account for the
!!!    Fermi energy
!!!
!!! Copyright 2009-2012 Philipp Wissgott
!!!           2014-2015 Elias Assmann

PROGRAM shift_energy
  use reallocate, only: realloc
  use clio,       only: argstr, fetcharg, croak
  use wien2k,     only: gtfnam, errflg, errclr
  use util,       only: string
  use const,      only: DPk, BUFSZ

  implicit none

  integer, parameter :: unit_eig=11, Ninit=10, Ninc=10

  real(DPk) :: EFermi

  real(DPk), allocatable :: energies(:)

  character(BUFSZ) :: deffn, errfn, status, form, fname

  integer :: kpt, band, num_bands, num_kpts
  integer :: iunit, ik, ib, iline, nbuf

  call fetcharg(2, EFermi)

  call gtfnam(deffn, errfn, ik)
  call errflg(errfn, 'Error in FINDBANDS')
  open (1, FILE=deffn, STATUS='old')
  do
     read(1, *, END=110) iunit, fname, status, form
     open(iunit, FILE=fname, STATUS=status, FORM=form)
  end do
110 continue

  allocate(energies(Ninit))
  nbuf = Ninit

  iline=1
  readeig: do
     read(unit_eig, *, END=120) kpt, band, energies(iline)
     if (iline >= nbuf) then
        call realloc(energies, nbuf+Ninc)
        nbuf = nbuf+Ninc
     end if
     iline = iline+1
  end do readeig
120 continue

  num_bands = band
  num_kpts  = kpt

  if (iline-1 /= num_bands*num_kpts) &
       call croak('num_bands='//trim(string(num_bands))// &
       &          ' and num_kpts='//trim(string(num_kpts))// &
       &          ' inconsistent with file size='//trim(string(iline-1)))

  rewind(unit_eig)

  writeeig: do ik = 1,num_kpts
     do ib = 1,num_bands
        write(unit_eig, '(2I12, F22.16)') &
             ib, ik, energies((ik-1)*num_bands + ib) + EFermi
     end do
  end do writeeig

  call errclr(errfn)
END PROGRAM shift_energy


!!/---
!! Local Variables:
!! mode: f90
!! End:
!!\---
!!
!! Time-stamp: <2015-05-23 19:58:48 elias>
