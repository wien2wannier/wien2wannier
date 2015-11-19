!!! wien2wannier/SRC_trig/find_bands.f
!!!
!!!    Finds the bands which are within an energy interval
!!!
!!! Copyright 2009-2012 Philipp Wissgott
!!!           2013-2014 Elias Assmann
!!!
!!! $Id: findbands.f 191 2014-03-03 09:16:46Z assmann $

program find_bands
  use const,  only: BUFSZ, DPk, Ryd_eV
  use wien2k, only: gtfnam, errflg, errclr

  implicit none

  integer, parameter :: unit_out1=51, unit_out=6, unit_in=5

  real(DPk) :: emin, emax, efermi, eband
  integer   :: ik=0, pos, dpos, iunit
  integer   :: nband, hiband, loband
  integer   :: himin=999, himax=0, lomin=999, lomax=0, nmin=999, nmax=0

  character(BUFSZ) :: deffn, errfn, fname, status, form
  ! FIXME: this is used to read ‘output1’
  ! if the lines of that file are too long, we fail silently
  character(BUFSZ) :: buf

  call gtfnam(deffn, errfn, ik)
  call errflg(errfn, 'Error in FINDBANDS')
  open (1, FILE=deffn, STATUS='old')
  do
     read(1, *, END=110) iunit, fname, status, form
     open(iunit, FILE=fname, STATUS=status, FORM=form)
  end do
110 continue 

  read(unit_in,*) emin, emax, efermi

  write(unit_out, '("# findbands: emin=", F8.3, " eV emax=", F8.3, &
       &" eV efermi=", F8.4, " Ry")') emin, emax, efermi
  write(unit_out, '(A)') '# k-point first last #bands'

  !read seed.output1
  kpoint: do
     ik = ik+1
     buf=''
     nband=0
     hiband=-1
     loband=1

     ! skip to eigenvalues
     do while (index(buf, 'EIGENVALUES ARE') == 0)
        read(unit_out1, '(A)', END=120) buf
     end do

     ! read eigenvalues
     line: do
        read(unit_out1, '(A)') buf
        if (index(buf, 'EIGENVALUES BELOW') /= 0) goto 130
        pos=1

        value: do
           dpos = scan(buf(pos:), '+-0123456789eE.')
           if (dpos==0) cycle line
           pos = pos + dpos-1

           read(buf(pos:), *) eband
           eband = (eband - efermi) * ryd_ev
           pos = pos + verify(buf(pos:), '+-0123456789eE.')

           if (eband > emax) goto 130

           if (eband >= emin) then 
              nband = nband+1
           else
              loband = loband+1
           end if
        end do value
     end do line

130  hiband=loband+nband-1
     write(unit_out, '(I6, "  ", I5, " ", I5, " ", I5)') &
          ik, loband, hiband, nband
     if (loband < lomin) lomin=loband
     if (loband > lomax) lomax=loband
     if (hiband < himin) himin=hiband
     if (hiband > himax) himax=hiband
     if (nband  < nmin)  nmin =nband
     if (nband  > nmax)  nmax =nband
  end do kpoint

120 continue

  write(unit_out,*)
  write(unit_out, '("Bloch bands in the interval")')
  write(unit_out, '("at all k: ", I3, " ", I5, " ", I5)') &
       lomax, himin, nmin
  write(unit_out, '("at any k: ", I3, " ", I5, " ", I5)') &
       lomin, himax, nmax

  call errclr(errfn)
end program find_bands

!!/---
!! Local Variables:
!! mode: f90
!! End:
!!\---
