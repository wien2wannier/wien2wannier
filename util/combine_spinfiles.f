PROGRAM combine_spinfiles
 !program combines amn,mmn files for up/dn spin from wien2k
 !to one amn,mmn by entrywise summation
 
 use util 

 implicit none
  character*50 seedname

  integer iarg,argcount !command line input argument counter
  integer :: nemin,nemax
  integer :: nproj  !#projection = # final Wannier functions
  integer :: nb !number of Bloch bands
  integer :: nk !number of kpoints
  integer :: nnnk !number of nearest neighbor kpoints
  integer :: idxk,idxnn !index variables for k points and nn k points
  integer :: idxb,idxp !index variables for bands and projections
  integer :: jk,jnn,jp,jb ! counter
  real*8 :: rvalup,rvaldn,cvalup,cvaldn !value storage
  integer :: bx,by,bz !storage for k vector
  character*20 title
  character*70 dummy
  
    
 !command line argument read-in
 iarg=iargc()
 argcount = 1
        
 if (argcount.eq.1) then
     call getarg(1,dummy)
     if (dummy(1:2).eq."-h") then
        write(*,*) 'program combines amnup/dn,mmnup/dn files from up/dn spin from wien2wannier'
        write(*,*) 'to one amn,mmn file via entrywise summation'
        write(*,*) 'Usage: combine_spinfiles case'
        stop
     else
        read(dummy,*)seedname
     endif
 else
     write(*,*) 'Usage: combine_spinfiles case'
     stop
 endif
 write(*,*)"case is:", seedname
 open(unit=1,file=clearspace(seedname)//'.mmnup',status='old')
 open(unit=2,file=clearspace(seedname)//'.mmndn',status='old')
 open(unit=3,file=clearspace(seedname)//'.mmnso',status='unknown')
 read(1,"(A20)") title
 read(1,"(3I12)")nb,nk,nnnk                                           
 read(2,"(A70)")dummy 
 read(2,"(A70)")dummy
 write(3,"(A20)") title
 write(3,"(3I12)")nb,nk,nnnk                                           

 do jk=1,nk

    do jnn=1,nnnk
            
          read(1,"(5i8)") idxk,idxnn,bx,by,bz
          read(2,"(A70)")dummy
          write(3,"(5i8)")idxk,idxnn,bx,by,bz
          do jb=1,nb**2
             read(1,"(32f18.12)")rvalup,cvalup
             read(2,"(32f18.12)")rvaldn,cvaldn
             write(3,"(32f18.12)")rvalup+rvaldn,cvalup+cvaldn
          enddo

    enddo

 enddo

 close(1)
 close(2)
 close(3)

 open(unit=1,file=clearspace(seedname)//'.amnup',status='old')
 open(unit=2,file=clearspace(seedname)//'.amndn',status='old')
 open(unit=3,file=clearspace(seedname)//'.amnso',status='unknown')
 read(1,"(A20)") title
 read(1,"(3I12)")nb,nk,nproj 
 read(2,"(A70)")dummy 
 read(2,"(A70)")dummy
 write(3,"(A20)") title
 write(3,"(3I12)")nb,nk,nproj

  do jk=1,nk
     do jp=1,nproj
	do jb=1,nb
	   read(1,"(2i4,1x,i5,1x,2e18.5)")idxb,idxp,idxk,rvalup,cvalup
           read(2,"(2i4,1x,i5,1x,2e18.5)")idxb,idxp,idxk,rvaldn,cvaldn
           write(3,"(2i4,1x,i5,1x,2e18.5)")idxb,idxp,idxk,rvalup+rvaldn,cvalup+cvaldn
        enddo
     enddo
  enddo

 close(1)
 close(2)
 close(3)
 write(*,*)"combined amnup+amndn and mmnup+mmndn files to amnso and mmnso"

 
END PROGRAM



