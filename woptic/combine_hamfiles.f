  PROGRAM combine_hamfiles
  !joins two Hamiltonian files and two mommat files for the woptic algorithm
  !the header of the first file is kept and the second on is appended to the
  !other.
  !It's case.ham_old + case.ham_fine = case.ham_joined
  !and  case.mommat_old + case.mommat = case.mommat_joined

  implicit none

  integer nk1,nk2,jk,nb1,nb2,nbmax,i
  integer hamunit1,hamunit2,iarg,j,ii,jj
  integer mommatunit1,mommatunit2,targethamunit,targetmomunit
  integer nemin,nemax
  real*8 qxte,qyte,qzte
  real*8, allocatable :: aa(:,:)
  complex*16 Xi
  parameter(nbmax = 99)
  parameter(Xi=(0.d0,1.d0))
  character*70 argdummy
  character*50 seedname
  character*120 dummy


  hamunit1 = 11
  hamunit2 = 12
  mommatunit1 = 13
  mommatunit2 = 14
  targethamunit = 15
  targetmomunit = 16
  iarg=iargc() 
   if(iarg.ge.1) then
      do j=1,iarg
         call getarg(j,argdummy)
         if (argdummy(1:1).eq.'-') then
            if ((argdummy(2:2).eq.'h').or.(argdummy(2:2).eq.'H')) then     
               write(*,*)"joins two Hamiltonian files and two mommat files for"
               write(*,*)"the woptic algorithm. The header of the first file is"
               write(*,*)"kept and the second on is appended to the other."
               write(*,*)"It's case.ham_old + case.ham_fine = case.ham_joined"
              write(*,*)"and case.mommat_old + case.mommat =case.mommat_joined"
               stop
            else
                 write(*,*)"Error: Unknown option"
                 stop
            endif
         else
            read(argdummy,*)seedname
         endif
      enddo

   else
      write(*,*) 'Usage: combine_hamfiles case'
      stop
   endif

  open(hamunit1,file=trim(seedname)//'.ham_old',status='old')
  open(hamunit2,file=trim(seedname)//'.ham_fine',status='old')
  open(targethamunit,file=trim(seedname)//'.ham_joined',status='unknown')
  read(hamunit1,*)nk1,nb1
  read(hamunit2,*)nk2,nb2
  if ((nb1.ne.nb2).and.(nk1.ne.0)) then
     write(*,*)"Error: number of bands not consistent"
     stop
  endif
  write(targethamunit,*)nk1+nk2,nb2
!   allocate(Hk(nk1,nb1,nb1),Hkdx(nk1,nbmax,nbmax))
  allocate(aa(nb2,2*nb2))
  do i=1,nk1
      read(hamunit1,*) qxte,qyte,qzte
      write(targethamunit,*) qxte,qyte,qzte
      do j=1, Nb2
         read(hamunit1,*) (aa(j,ii),ii=1,2*Nb2)
      enddo
      do j=1, Nb2
         write(targethamunit,*) (aa(j,ii),ii=1,2*Nb2)
      enddo
  enddo
  do i=1,nk2
      read(hamunit2,*) qxte,qyte,qzte
      write(targethamunit,*) qxte,qyte,qzte
      do j=1, Nb2
         read(hamunit2,*) (aa(j,ii),ii=1,2*Nb2)
      enddo
      do j=1, Nb2
         write(targethamunit,*) (aa(j,ii),ii=1,2*Nb2)
      enddo
  enddo 
  close(hamunit1)
  close(hamunit2)
  close(targethamunit)

  open(mommatunit1,file=trim(seedname)//'.mommat_old',status='old')
  open(mommatunit2,file=trim(seedname)//'.mommat',status='old')
  open(targetmomunit,file=trim(seedname)//'.mommat_joined',status='unknown')
  read(mommatunit1,"(A120)")dummy
  write(targetmomunit,"(A120)")dummy
  do j=1,nk1
      read(mommatunit1,"(A120)")dummy
      write(targetmomunit,"(A120)")dummy
      read(mommatunit1,"(A120)")dummy
      write(targetmomunit,"(A120)")dummy
      read(dummy(28:38),"(2I5)")nemin,nemax
      read(mommatunit1,"(A120)")dummy
      write(targetmomunit,"(A120)")dummy
      do ii=nemin,nemax
           do jj=nemin,nemax
               if (jj.ge.ii) then
                     read(mommatunit1,"(A120)")dummy
                     write(targetmomunit,"(A120)")dummy
!                      read(dummy,9040)i1,i2,mxr,mxi,myr,myi,mzr,mzi,ediff
               endif
           enddo
      enddo
  enddo
  read(mommatunit2,"(A120)")dummy
  do j=1,nk2
      read(mommatunit2,"(A120)")dummy
      write(targetmomunit,"(A120)")dummy
      read(mommatunit2,"(A120)")dummy
      write(targetmomunit,"(A120)")dummy
      read(dummy(28:38),"(2I5)")nemin,nemax
      read(mommatunit2,"(A120)")dummy
      write(targetmomunit,"(A120)")dummy
      do ii=nemin,nemax
           do jj=nemin,nemax
               if (jj.ge.ii) then
                     read(mommatunit2,"(A120)")dummy
                     write(targetmomunit,"(A120)")dummy
!                      read(dummy,9040)i1,i2,mxr,mxi,myr,myi,mzr,mzi,ediff
               endif
           enddo
      enddo
  enddo
  close(mommatunit1)
  close(mommatunit2)
  close(targetmomunit)



  END PROGRAM