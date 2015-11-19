  PROGRAM combine_plots

  use util

  implicit none

  integer iarg,j,jb
  character*80 argdummy,file1,file2,file3,file4
  character*180 dummy
  integer np(3)
  real*8 :: tmpr(4,10),resv(10),resp(10)
  complex*16 :: tmpc,Xi

  Xi = dcmplx(0d0,1d0)
  iarg=iargc()
   if ((iarg.ge.1)) then
    do j=1,iarg
        call getarg(j,argdummy)
        if (argdummy(1:1).eq.'-') then
!             if ((argdummy(2:3).eq.'up').or.(argdummy(2:3).eq.'dn')) then     
!                !for spin-polarized calc. the fileendings have additional up/dn
!               vectorfileend = ".vector"//argdummy(2:3)
!               energyfileend = ".energy"//argdummy(2:3)
!               startmessage = "++ join vector files of spin-polarized input files:"//argdummy(2:3)//" ++"
!             elseif ((argdummy(2:5).eq.'soup').or.(argdummy(2:5).eq.'sodn')) then  
!                vectorfileend = ".vectorso"//argdummy(4:5)
!                energyfileend = ".energyso"//argdummy(4:5)
!                usecomplex = .true.
!                startmessage = "++ join vector files of spin-polarized spin-orbit input files:"//argdummy(4:5)//" ++"
!             elseif (argdummy(2:2).eq.'c') then  
!                usecomplex = .true.
             if (argdummy(2:2).eq.'h') then
                write(*,*) 'joins two plot files plot1 and plot2 (files should be psink files)'
                write(*,*) 'by adding them considering the corresponding phases in the psiarg files'
                write(*,*) 'Usage: combine_plots file1(psink file) file2(psiarg file from file1)'// & 
                            'file3(psink file) file4(psiarg file from file3)'
                write(*,*) 'the names of the outputfiles are combined.psink and combined.psiarg'
                stop    
             else
                write(*,*) 'Unknown option1'
                write(*,*) 'Usage: combine_plots file1(psink file) file2(psiarg file) file3(psink file) file4(psiarg file)'
                stop    
             endif
         else
            if ((iarg.ne.4)) then
                write(*,*) 'Unknown option1'
                write(*,*) 'Usage: combine_plots file1(psink file) file2(psiarg file) file3(psink file) file4(psiarg file)'
                stop    
            endif
            if (j.eq.1) then
               read(argdummy,*)file1
            elseif (j.eq.2) then
               read(argdummy,*)file2
            elseif (j.eq.3) then
               read(argdummy,*)file3
            elseif (j.eq.4) then
               read(argdummy,*)file4   
            endif
        endif
     enddo
 else
   write(*,*) 'xxUsage: combine_plots file1(psink file) file2(psiarg file) file3(psink file) file4(psiarg file)'
   stop
 endif


  open(unit=1,file=clearspace(file1),status='old')
  open(unit=2,file=clearspace(file2),status='old')
  open(unit=3,file=clearspace(file3),status='old')
  open(unit=4,file=clearspace(file4),status='old')
  open(unit=10,file='combined.psink',status='unknown')
  open(unit=11,file='combined.psiarg',status='unknown')
  read(1,"(A180)")dummy
  write(10,"(A180)")dummy
  write(*,*)dummy
  read(1,"(A180)")dummy
  write(10,"(A180)")dummy
  read(dummy,*)np(1)
  read(1,"(A180)")dummy
  write(10,"(A180)")dummy
  read(dummy,*)np(2)
  read(1,"(A180)")dummy
  write(10,"(A180)")dummy
  read(dummy,*)np(3)
  read(1,"(A180)")dummy
  write(10,"(A180)")dummy
  do j=1,5
     read(3,*)dummy
  enddo
  do j=1,np(1)*np(2)*np(3)/10
     !readin
     read(1,*)tmpr(1,:)
     read(2,*)tmpr(2,:)
     read(3,*)tmpr(3,:)
     read(4,*)tmpr(4,:)
     !combine
     do jb=1,10
        tmpc = dcmplx(dsqrt(tmpr(1,jb)))*cdexp(Xi*tmpr(2,jb)) + &
               dcmplx(dsqrt(tmpr(3,jb)))*cdexp(Xi*tmpr(4,jb))
        resv(jb) = cdabs(tmpc)**2
        resp(jb) = datan2(dimag(tmpc),dble(tmpc))
     enddo
     write(10,3020)resv
     write(11,3020)resp
  enddo
  close(1)
  close(2)
  close(3)
  close(4)
  close(10)
  close(11)
  write(*,*)file1,np
  write(*,*)file2
  write(*,*)file3
  write(*,*)file4


3020 FORMAT(1P,10E16.8)
  END PROGRAM