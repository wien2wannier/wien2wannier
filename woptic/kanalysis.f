program kanalysis

use util

implicit none

character*80 seedname
character*80 argdummy,dumrest
integer nkloc,nk,jk,iarg,unit1,unit2
integer niter,jf,nw,dumi,kdiv,kloc(3),ik,dumi3
integer jw,nkband,j1,counter,jk1,jk2,j2,js,kmaxidx
integer, allocatable :: kminidx(:)
real*8, allocatable :: kband(:,:),kp(:,:),kwdata(:,:,:)
real*8 :: kdist,kdistloc,kdistsum,rvec(3),kmax,convfac,dumr,dw
character*1 dum1,c1
character*2 c2
character*10 dum10
type(structure) lattice
logical matlabmode

  iarg=iargc() 
  if (iarg.lt.2) then
     write(*,*)"Usage: kanalysis startingfrequency case mode(optional)"
     stop
  endif
  matlabmode = .false.
   if (iarg.eq.3) then
       call getarg(3,argdummy)
       if (argdummy(1:1).eq."1") then
          matlabmode = .true.
          write(*,*)"Matlab mode"
       endif
   endif
  
  call getarg(2,argdummy)
  write(seedname,"(A80)")argdummy
  call getarg(1,argdummy)
  read(argdummy,*)niter



  unit1 = 1
  open(unit=unit1,file=clearspace(seedname)//'.woptin',status='old')
  read(unit1,*)
  read(unit1,*)dumr,dw
  close(unit1)
  unit2 =2
  open(unit=unit1,file=clearspace(seedname)//'.struct',status='old')
  open(unit=unit2,file=clearspace(seedname)//'.outputkgen_orig',status='old')
  call countatoms(unit1,lattice)
  call readin_lattice(unit1,unit2,lattice)
   !round transform matrix
     do j1=1,3
     do j2=1,3
        lattice%transform_matrix(j1,j2) = round_real(lattice%transform_matrix(j1,j2),4)
     enddo
  enddo
   close(unit1)
   close(unit2)

   open(unit=unit1,file=clearspace(seedname)//'.kcontribw_band',status='old')
   open(unit=unit2,file=clearspace(seedname)//'.optanalysis_band',status='unknown')
   read(unit1,*)dum1,nk,nw,dumi,convfac
   write(unit2,*)dum1,nk,nw+1-niter,convfac
   kmax = 0d0
   allocate(kwdata(nk,0:nw,6),kp(nk,3))
    do jk=1,nk
         read(unit1,*)dumi,kwdata(jk,0,:)
         do jw=1,nw
             read(unit1,*)kwdata(jk,jw,:)
             if (jw.ge.niter) then
                 write(unit2,*)jk,dw*jw,kwdata(jk,jw,1)*convfac/1000d0
             endif
         enddo
         if (kwdata(jk,150,1).gt.kmax) then
               kmax = kwdata(jk,150,1)
               kmaxidx = jk
        endif
        if  (.not.matlabmode) then
           write(unit2,*)
        endif
    enddo
    write(*,*)'maximumk',kmax,kmaxidx
!   stop
   close(unit1)
   open(unit=unit1,file=clearspace(seedname)//'.K1w_band',status='old')
   open(unit=unit2,file=clearspace(seedname)//'.thermanalysis_band',status='unknown')
!    read(unit1,*)dum1,nk,nw,dumi,convfac
   write(unit2,*)dum1,nk,nw+1,convfac
   kmax = 0d0
    do jk=1,nk
!          read(unit1,*)dumi,kwdata(jk,0,:)
!          do jw=1,nw
             read(unit1,*)dumi,kwdata(jk,0,:)
              write(unit2,*)jk,kwdata(jk,0,1)
!          enddo
         if (kwdata(jk,0,1).gt.kmax) then
               kmax = kwdata(jk,0,1)
               kmaxidx = jk
        endif
!         if  (.not.matlabmode) then
!            write(unit2,*)
!         endif
    enddo
    write(*,*)'maximumk',kmax,kmaxidx
!   stop
   close(unit1)


end program