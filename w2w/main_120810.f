      PROGRAM WF                                                      
! LAPW2DM calculates the density matrix
! for new calculation of <X> and for rotationally invariant LDA+U method
! LAPW2DM is modified LAPW2 packages
! last change P.Novak 08.010.2001 novak fzu.cz
        USE param
        USE struct
        USE xa
        USE xa3
        USE bessel
        USE amn_mmn
        USE ams
        USE pairs
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16      CIM
      CHARACTER*10    KNAME
      CHARACTER*3     MODE
      CHARACTER*11     STATUS,FORM                                      
      CHARACTER*67       ERRMSG
      CHARACTER*80       DEFFN, ERRFN,aline
      CHARACTER*180     FNAME,VECFN  
      LOGICAL   MMN,AMN
      INTEGER CENTERATOM(100)

      COMMON /GENER/   BR1(3,3),BR2(3,3)     
      DATA  CIM /(0.d0,1.d0)/
                            
!                                                                       
      
!-----------------------------------------------------------------------  
!                                                                       
      CALL init_ams
      CALL GTFNAM(DEFFN,ERRFN,IPROC)
      CALL ERRFLG(ERRFN,'Error in LAPW2DM')
      OPEN (1,FILE=DEFFN,STATUS='OLD',ERR=910)
   10 CONTINUE
         READ (1,*,END=20,ERR=960) IUNIT,FNAME,STATUS,FORM,IRECL
         OPEN (IUNIT,FILE=FNAME,STATUS=STATUS,FORM=FORM,ERR=920)
      GOTO 10
   20 CONTINUE
      CLOSE (1)

!.....READ STRUCT                                                       
      CALL init_struct
!....Find nume, nmat and nkpt
     k=0  
      DO I=1,NAT 
         READ(50,'(f9.5)') EMIST
         READ(50,'(f9.5)') EMIST
      ENDDO
!      DO I=1,NAT                                                  
!         READ(51,'(f9.5)',iostat=ios) EMIST
!         READ(51,'(f9.5)',iostat=ios) EMIST
!      ENDDO
!      IF(ios.eq.0) JSPIN=2
      DO
         READ(50,'(3e19.12,a10,2i6)',IOSTAT=ios) Sxx,Txx,Zxx,KNAME,N,NEn
         IF (ios /= 0) EXIT  
         k=k+1
         nmat=MAX(n,nmat)
         nume=MAX(nen,nume)
         DO ii=1,nen
            READ(50,*) NUM,E1
         ENDDO
!         IF(jspin.EQ.2) THEN
!            READ(51,'(3e19.12,a10,2i6)',IOSTAT=ios) Sxx,Txx,Zxx,KNAME,N,NEn
!            nmat=MAX(n,nmat)
!            nume=MAX(nen,nume)
!            DO ii=1,nen
!               READ(51,*) NUM,E1
!            ENDDO
!         ENDIF
      ENDDO
      nkpt=k+1
      REWIND(50)
!      REWIND(51)
!.....READ INPUT AND POTE  
! if more information needed, put IPRINT>0
      IPRINT=0
 113  read(10,'(a80)',end=114) aline
      if(index(aline(1:80),'begin kpoints').ne.0) then
	read(10,*,err=913)NUM_KPTS 
      endif
      if(index(aline(1:80),'begin nnkpts').ne.0) then
	read(10,*,err=914)NNTOT
        N_pair=NUM_KPTS*NNTOT
        CALL init_pairs(n_pair)
!        if (N_pair.gt.NPAIR) stop 'NPAIR too small check modules'
        do i=1,N_pair
         READ(10,*)KP(I),KPB(I),BQX1(I),BQY1(I),BQZ1(I)
        enddo
       endif
       goto 113
 114   continue
       write(6,*)'NUM_KPTS=',NUM_KPTS
       write(6,*)'NNTOT=',NNTOT
       write(6,*)'N_pair=',N_pair

      MMN=.true.
      AMN=.true.
      READ(5,*)MODE
      IF (MODE.eq.'MMN') AMN=.false.  
      IF (MODE.eq.'AMN') MMN=.false.
      write(6,*)'MODE=' 
      if (MMN) write(6,*)'MMN'
      if (AMN) write(6,*)'AMN'
      READ(5,*)NEMIN,NEMAX
	if (nemax.lt.nemin) stop 'nemin > nemax'
      READ(5,*)LJMAX,NPROJ
      write(6,*)'nemin,nemax:',nemin,nemax
      nb=nemax-nemin+1
      CALL init_xa3(nb,nmat,nkpt)
      CALL init_xa(LMAX2,NMAT,NRAD,nb)
      CALL read_vec(nemin,nemax,nkpt,maxx,maxy,maxz)
      CALL init_mmn(nb,n_pair,nproj)
       
      C(:,:,:) = cmplx(0d0,0d0)

      IF (AMN) THEN
       DO I=1,NPROJ
        READ(5,*,err=925,end=926)N
        DO J=1,N
         READ(5,*,err=925,end=926)IA,L,M,X1,X2
         ind=L*(L+1)+M+1
         C(I,ind,IA)=X1+CIM*X2
         CENTERATOM(I)=IA 
        ENDDO
       ENDDO
       write(6,*)'initial orbital projections'
       DO I=1,nproj
        write(6,*)'orbital #',I,'centered at atom',CENTERATOM(I)
        ind=0
        DO L=0,3
         DO M=-L,L
          ind=ind+1
          write(6,*)L,M,C(I,ind,CENTERATOM(I))
         ENDDO
        ENDDO
       ENDDO
      ENDIF

      call init_bessel(LMAX2,LJMAX,NRAD,NRF)
      call gaunt2
 1234 FORMAT(//,1A)
      WRITE(6,800)                                                      
      WRITE(6,805)  TITLE    
      if (amn) then
       write(7,806) TITLE 
       write(7,807) NEMAX-NEMIN+1,NUM_KPTS,NPROJ
      endif
      if (mmn) then
       write(8,806) TITLE
       write(8,807) NEMAX-NEMIN+1,NUM_KPTS,NNTOT                                           
      endif
      WRITE(6,810)  LATTIC                                              
      WRITE(6,820)  AA,BB,CC                                            
      WRITE(6,840)  NAT                                                 
      WRITE(6,850)  IREL                                                
!      WRITE(6,870)  COORD  
      CALL LATGEN
!     rotate boundary vectors
       DO I=1,N_pair
        BX=BQX1(I)
        BY=BQY1(I)
        BZ=BQZ1(I)
        BQX(I)=BX*BR2(1,1)+BY*BR2(1,2)+BZ*BR2(1,3)    
        BQY(I)=BX*BR2(2,1)+BY*BR2(2,2)+BZ*BR2(2,3)
       BQZ(I)=BX*BR2(3,1)+BY*BR2(3,2)+BZ*BR2(3,3)
       ENDDO
      write(6,*)' alpha test',(alpha(i),i=1,3)
!.....CALCULATE CHARGE DENSITY CLM(R) IN SPHERES,  PARTIAL CHARGES      

      if (MMN) then
       call cputim(t1)
       call l2mmn(NB,NUM_KPTS,NNTOT,LJMAX)
       call cputim(t2)
       !JXZ: MAXG is not pre-defined
        MAXG = 0
        write(6,*)'MXG=',MAXG
       call planew(NB,NUM_KPTS,NNTOT,maxx,maxy,maxz)
       call cputim(t3)
       write(6,*)'CPU l2mmn:',t2-t1
       write(6,*)'CPU planew:',t3-t2
      endif
      if (AMN) then
       call cputim(t1)
       call l2amn(NB,NPROJ,NUM_KPTS)
       call cputim(t2)
       write(6,*)'CPU l2amn:',t2-t1
      endif

      CALL ERRCLR(ERRFN)
      STOP 'W2W END'                                                 
!
!        error handling
!
  910 INFO = 1
!
!        'lapw2.def' couldn't be opened
!
      WRITE (ERRMSG,9000) FNAME
      CALL OUTERR('lapwdm',ERRMSG)
      GOTO 999
  920 INFO = 2
!
!        file FNAME couldn't be opened
!
      WRITE (ERRMSG,9010) IUNIT
      CALL OUTERR('lapwdm',ERRMSG)
      WRITE (ERRMSG,9020) FNAME
      CALL OUTERR('lapwdm',ERRMSG)
      WRITE (ERRMSG,9030) STATUS, FORM
      CALL OUTERR('lapwdm',ERRMSG)
      GOTO 999
  925 WRITE(*,*)'error reading projection',I
      GOTO 999
  926 WRITE(*,*)'too few projections'
  930 INFO = 3
!
!        illegal number of equivalent atoms
!
      CALL OUTERR('lapwdm','MULT .EQ. 0')
      GOTO 999
  950 INFO = 5
      CALL OUTERR('lapwdm','Too many atoms (NATO too small).')
      GOTO 999
  955 INFO = 6
      CALL OUTERR('lapwdm','Too many atoms (NDIF too small).')
      GOTO 999
  956 INFO = 56
      CALL OUTERR('lapwdm','LXDOS must be 3 for ISPLIT=999.')
      GOTO 999
  960 INFO = 7
  913 CALL OUTERR('error reading begin kpoints')
      GOTO 999
  914 CALL OUTERR('error reading begin nnkpts')
      GOTO 999
!
!        Error reading file 'lapw2.def'
!
      WRITE (ERRMSG,9040) FNAME
      CALL OUTERR('lapwdm',ERRMSG)
      GOTO 999
  999 STOP 'lapwdm - Error'
!                                                                       
!                                                                       
  43  FORMAT(3X,A77)                                                    
  44  FORMAT(I3,A77)                                                    
 700  FORMAT(I3,A77)                                                    
 800  FORMAT(////,30X,50(1H-),/,33X,'S T R U C T U R A L   ',            &
             'I N F O R M A T I O N',/,30X,50(1H-),//)                  
 805  FORMAT(3X,'SUBSTANCE',20X,'= ',A80,/)   
 806  FORMAT(A20)  
 807  FORMAT(3I12)                        
 810  FORMAT(3X,'LATTICE',22X,'= ',A4)                                  
 820  FORMAT(3X,'LATTICE CONSTANTS ARE',8X,'= ',3F12.7)                 
 830  FORMAT(3X,'SYMMETRY ATTITUDE IS',9X,'= ',A4)                      
 840  FORMAT(3X,'NUMBER OF ATOMS IN UNITCELL  = ',I3)                   
 850  FORMAT(3X,'MODE OF CALCULATION IS',7X,'= ',A4)                    
 860  FORMAT(3X,'SELFCONSISTENT CYCLE-NUMBER  = ',I3,/)                 
 870  FORMAT(3X,'TYPE OF COORDINATES IN DSPLIT= ',A5)                   
 1000 FORMAT(A80)                                                       
 1010 FORMAT(A4,24X,I2,1x,a4,/,13X,A4,18X,A4)                                 
 1020 FORMAT(6F10.7,10X,F10.7)                                          
 1030 FORMAT(5X,I3,4X,F10.7,3X,F10.7,3X,F10.7,/,15X,I2,17X,I2)          
 1031 FORMAT(5X,I3,4X,F10.7,3X,F10.7,3X,F10.7)                          
 1040 FORMAT(///,3X,'ERROR IN lapwdm : MULT(JATOM)=0 ...',                &
             /, 20X,'JATOM=',I3,3X,'INDEX=',I3,3X,'MULT=',I3)           
 1050 FORMAT(A10,5X,I5,5X,F10.9,5X,F10.5,5X,F10.5)                      
 1002 FORMAT(3F10.5,I5)                                                 
 1003 FORMAT(A5)                                                        
 1004 FORMAT(A5,f10.5)                                                        
 1060 FORMAT(//,':NOE  :',1X,'NUMBER OF ELECTRONS',10X,'= ',F7.3,            &
             //   ,':FER  :',1X,'F E R M I - ENERGY',11X,'= ',F9.5)            
 1061 FORMAT(//,':NOE  :',1X,'NUMBER OF ELECTRONS',10X,'= ',F7.3,            &
       //   ,':FER  :',1X,'F E R M I - ENERGY(TETRAH.M.)','= ',F9.5)           
 1062 FORMAT(//,':NOE  :',1X,'NUMBER OF ELECTRONS',10X,'= ',F7.3,            &
       //   ,':FER  :',1X,'F E R M I - ENERGY(GAUSS-.M.)','= ',F9.5)           
 1063 FORMAT(//,':NOE  :',1X,'NUMBER OF ELECTRONS',10X,'= ',F7.3,            &
       //   ,':FER  :',1X,'F E R M I - ENERGY(FERMI-SM.)','= ',F9.5)           
 2000 FORMAT(//,3X,'=====>>> CPU TIME SUMMARY',/)                       
 2010 FORMAT(12X,'TOTAL       : ',F8.1,5X,'... ',F4.0,' PERCENT')       
 2020 FORMAT(12X,'PART FERMI  : ',F8.1,5X,'... ',F4.0,' PERCENT')       
 2030 FORMAT(12X,'PART CLM    : ',F8.1,5X,'... ',F4.0,' PERCENT')       
 2040 FORMAT(12X,'PART FOURIR : ',F8.1,5X,'... ',F4.0,' PERCENT')       
 6000 FORMAT(///,3X,'ERROR IN lapwdm : MULT(JATOM)=0 ...', &
             /, 20X,'JATOM=',I3,3X,'INDEX=',I3,3X,'MULT=',I3)
 9000 FORMAT('can''t open definition file ',A40)
 9010 FORMAT('can''t open unit: ',I2)
 9020 FORMAT('       filename: ',A50)
 9030 FORMAT('         status: ',A,'  form: ',A)
 9040 FORMAT('Error reading file: ',A47)
      END                                                               
