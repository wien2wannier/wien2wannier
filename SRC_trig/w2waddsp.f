!!! wien2wannier/SRC_trig/w2waddsp.f
!!!
!!!    Combines amn,mmn files for up/dn spin from wien2k to one
!!!    amn,mmn by entrywise summation
!!!
!!! Copyright 2009-2012 Philipp Wissgott
!!!           2014-2015 Elias Assmann

PROGRAM combine_spinfiles
  use const,  only: DPk, bufsz
  use wien2k, only: gtfnam, errflg, errclr
  use util,   only: string

  implicit none

  integer, parameter :: unit_mmn=13, unit_mmn1=11, unit_mmn2=12
  integer, parameter :: unit_amn=23, unit_amn1=21, unit_amn2=22

  integer          :: num_bands, num_kpts, nntot, num_wann
  real(DPk)        :: val1(2), val2(2)
  integer          :: idx1(3), idx2(3)
  character(BUFSZ) :: buf, deffn, errfn, status, form, fname

  integer :: ik, iunit, iline, iblock
  logical :: doMmn, doAmn

  call gtfnam(deffn, errfn, ik)
  call errflg(errfn, 'Error in FINDBANDS')
  open (1, FILE=deffn, STATUS='old')
  do
     read(1, *, END=110) iunit, fname, status, form
     open(iunit, FILE=fname, STATUS=status, FORM=form)
  end do
110 continue

  inquire(unit_mmn1, OPENED=doMmn)
  inquire(unit_amn1, OPENED=doAmn)

!!! case.mmn
  if (doMmn) then
     read (unit_mmn1, *)         buf
     read (unit_mmn2, *)
     write(unit_mmn, '(A)') trim(buf)
     read (unit_mmn1, *)        num_bands, num_kpts, nntot
     read (unit_mmn2, *)
     write(unit_mmn,  '(3I12)') num_bands, num_kpts, nntot
     mmnblock: do iblock = 1, num_kpts*nntot
        read (unit_mmn1, '(A)')     buf
        read (unit_mmn2, *)
        write(unit_mmn, '(A)') trim(buf)

        mmnline: do iline = 1,num_bands**2
           read (unit_mmn1, *)          val1
           read (unit_mmn2, *)          val2
           write(unit_mmn, '(2F18.12)') val1+val2
        end do mmnline
     end do mmnblock
  end if

!!! case.amn
  if (doAmn) then
     read (unit_amn1, *)         buf
     read (unit_amn2, *)
     write(unit_amn, '(A)') trim(buf)
     read (unit_amn1, *)        num_bands, num_kpts, num_wann
     read (unit_amn2, *)
     write(unit_amn,  '(3I12)') num_bands, num_kpts, num_wann
     amnline: do iline = 1, num_bands*num_kpts*num_wann
        read (unit_amn1, *)                      idx1, val1
        read (unit_amn2, *)                      idx2, val2
        write(unit_amn, '(2I4,1X,I5,1X,2E18.5)') idx1, val1+val2
     end do amnline
  end if

  call errclr(errfn)
END PROGRAM combine_spinfiles


!!/---
!! Local Variables:
!! mode: f90
!! End:
!!\---
!!
