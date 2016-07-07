!!! wien2wannier/SRC_w2w/mknam.f

subroutine MKNAM(FNAME,OLDNAM,ILOOP)
! create a filename with running index ILOOP from
! a given 'parental' file name OLDNAM

  implicit none

  character(180) :: FNAME, OLDNAM
  character(4)   :: ALOOP

  integer :: iloop, ifrom, i, ito

  WRITE(ALOOP,'(I4)')ILOOP

  FNAME=OLDNAM

  IFROM=1
  do I=len(ALOOP),1,-1
     if(ALOOP(I:I) /= ' ') IFROM=I
  enddo

  ITO=len(ALOOP)

  do I=len(FNAME),1,-1
     if(FNAME(I:I).ne.' ') then
        FNAME(I+1:len(FNAME))='_' // ALOOP(IFROM:ITO)
        goto 7778
     endif
  end do
7778 continue
end subroutine MKNAM


!!/---
!! Local Variables:
!! mode: f90
!! End:
!!\---
!!
!! Time-stamp: <2016-07-07 12:26:26 assman@faepop71.tu-graz.ac.at>
