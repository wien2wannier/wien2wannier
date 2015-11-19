!$hp9000_800 intrinsics on
      PROGRAM LAPW7
!     last changes: 29.08.00 ub (removing CALL to PRTRH)
!
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*80      FNAME
      CHARACTER*11      STATUS,FORM
      
      !argument read-in parameters
      integer :: iarg,idx_wann
      character*3 dummy
!
!     << establish file conections >>
      iarg=iargc()
      CALL GETARG(1,FNAME)
      

      !command line optional 2nd argument read-in
       iarg=iargc()
       !argcount = 1
       if(iarg.eq.2) then
          call getarg(2,dummy)
          read(dummy,*)idx_wann
       else 
          idx_wann = 0
       endif
!      CALL GETARG(2,FNAME)
!      IF(FNAME.EQ.'      ') CALL GETARG(1,FNAME)
      OPEN(1,FILE=FNAME,STATUS='OLD',ERR=8000)
 8003 READ(1,*,END=8001) IUNIT,FNAME,STATUS,FORM,IRECL
      OPEN(IUNIT,FILE=FNAME,STATUS=STATUS,FORM=FORM,ERR=8002)
      GOTO 8003
!
 8000 WRITE(*,*) ' ERROR IN OPENING WPLOT.DEF !!!!'
      STOP 'WPLOT.DEF'
!
 8002 WRITE(*,*) ' ERROR IN OPENING UNIT:',IUNIT
      WRITE(*,*) '       FILENAME: ',FNAME,'  STATUS: ',STATUS, &
                 '  FORM:',FORM
      STOP 'OPEN FAILED'
!
 8001 CONTINUE
!
!     << start wave function evaluation >>
      CALL MAIN2(idx_wann)
      STOP
      END
