      PROGRAM kconv                                                      
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*80 aline,DEFFN,ERRFN,NNKP,KLIST

      CALL GTFNAM(DEFFN,ERRFN,IPROC)
      write(*,*)DEFFN
      NNKP=DEFFN
      KLIST=DEFFN
      DO 7777 I=LEN(DEFFN),1,-1
         IF(DEFFN(I:I).NE.' ') THEN
            KLIST(I+1:LEN(DEFFN))='.klist'  
            NNKP(I+1:LEN(DEFFN))='.nnkp'
            GOTO 7778
         ENDIF
 7777 CONTINUE
 7778 CONTINUE
      write(*,*)KLIST,NNKP
      OPEN (10,FILE=NNKP,STATUS='OLD',ERR=910)
      OPEN (11,FILE=KLIST,STATUS='UNKNOWN',ERR=920)



      write(*,*)'KFAC,emin,emax?'
      read(*,*)KFAC,emin,emax
 113  read(10,'(a80)',end=114) aline
      if(index(aline(1:80),'begin kpoints').ne.0) then
	read(10,*,err=913)NUM_KPTS 
        do K=1,NUM_KPTS
         read(10,*)x,y,z
	 kx=int(x*kfac)
         ky=int(y*kfac)
         kz=int(z*kfac)
	 if (k.eq.1) then
	  write(11,5)k,kx,ky,kz,kfac,dum,emin,emax
	else
	  write(11,6)k,kx,ky,kz,kfac,dum
	endif
       enddo
	write(11,7)'END'
      endif
      goto 113
 114  stop 'k-file converted'
 913  write(*,*)'error readin num_kpts'
      goto 999
 910  write(*,*)'error opening nnkp'
      goto 999
 920  write(*,*)'error opening klist'
      goto 999
     
   5  format(4i10,i9,3f5.1)
   6  format(4i10,i9,f5.1)
   7  format(a3)
 999     end

      SUBROUTINE GTFNAM(DEFFN,ERRFN,IPROC)
      CHARACTER*(*)      DEFFN, ERRFN
      INTEGER            IPROC
!        Local Parameters
!
      CHARACTER*5        ERREXT
      PARAMETER          (ERREXT = 'error')
!
!        Local Scalars
!
      INTEGER            I
!
!        extract the command-line argument
!
      IPROC=0
      iarg=iargc()
      if(iarg.eq.1) then
      CALL GETARG(iarg,DEFFN)
      else if(iarg.eq.2) then
         CALL GETARG(2,DEFFN)
            READ(DEFFN,*)IPROC
         CALL GETARG(1,DEFFN)
      else
  900 STOP 'GTFNAM - One or two commandline arguments have to be given.'
      endif
!
!        generate a name for the error-message file
!
      DO 10 I = LEN(DEFFN), 1, -1
         IF (DEFFN(I:I) .EQ. '.') THEN
            IF (LEN(ERRFN) .LT. (I+LEN(ERREXT))) GOTO 910
            ERRFN(1:I) = DEFFN(1:I)
            ERRFN(I+1:LEN(ERRFN)) = ERREXT
            GOTO 30
         ENDIF
   10 CONTINUE
!
!        the name of the definition file contains no '.', it is assumed
!        that this name contains no extension - append the extension
!        '.error' to get a name for the error file.
!
      DO 20 I = LEN(DEFFN), 1, -1
         IF (DEFFN(I:I) .NE. ' ') THEN
            IF (LEN(ERRFN) .LT. (I+1+LEN(ERREXT))) GOTO 910
            ERRFN(1:I) = DEFFN(1:I)
            ERRFN(I+1:LEN(ERRFN)) = '.' // ERREXT
            GOTO 30
         ENDIF
   20 CONTINUE
!
!        filename contains only spaces
!
      STOP 'GTFNAM - string ERRFN contains just spaces.'
   30 CONTINUE
!
      RETURN
      RETURN
!
!        Errors
!
  910 STOP 'GTFNAM - string ERRFN too short to hold filename.'
 920  STOP 'GTFNAM - number of parallel processes erroneous.'
!
!        End of 'GTFNAM'
!
      END


