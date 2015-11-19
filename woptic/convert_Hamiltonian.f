  PROGRAM convert_Hamiltonian
  !program to Fourier transform the real space Hamiltonian given in case_hr.dat
  !to k space H(k) where the k list is given by case.klist

  use util

  implicit none

  integer jk,jr,nb,nr,jb1,jb2,nk,nkp,kdim(3)
  integer hamunit,iarg,j,i,jd1,jd2
  integer targethamunit,klistunit,unitstruct,unitkgen
  integer di1,di2,IK!,line_count
  integer,allocatable :: rweights(:),rvec(:,:),k(:,:),kdiv(:)
  real*8 realval,imagval,emin,emax,tmpvec(3),pi
  real*8, allocatable :: wei(:),k2(:,:),loc_ham(:)
  complex*16 Xi,fac,rdotk
!   parameter(nbmax = 99)
  parameter(Xi=(0.d0,1.d0))
   parameter (Pi = dacos(-1.d0))
  character*70 argdummy
  character*50 seedname
  character*120 dummy
  complex*16, allocatable :: Hr(:,:,:),Hk(:,:,:)
  type(structure) lattice
  logical :: lattconv
  character*21 header
  character*10 dum10
  character*80 dumrest

  lattconv = .false.
  hamunit = 11
  targethamunit = 12
  klistunit = 13
  unitstruct = 14
  unitkgen = 15
  iarg=iargc() 
   if(iarg.ge.1) then
      do j=1,iarg
         call getarg(j,argdummy)
         if (argdummy(1:1).eq.'-') then
            if ((argdummy(2:2).eq.'h').or.(argdummy(2:2).eq.'H')) then     
               write(*,*)"program to Fourier transform the real space"         
               write(*,*)"Hamiltonian given in case_hr.dat to k space "
               write(*,*)"H(k) where the k list is given by case.klist"
               write(*,*)"flag -l leads to klist treatment as in H,R"
               write(*,*)"lattices"
               stop
            elseif (argdummy(2:2).eq.'l') then
               lattconv = .true.
            else
               write(*,*)"Unknown option"
               stop
            endif
         else
            read(argdummy,*)seedname
         endif
      enddo

   else
      write(*,*) 'Usage: convert_Hamiltonian [-l] case'
      stop
   endif

  !read-in from seed.struct/seed.outputkgen
   open(unit=unitstruct,file=clearspace(seedname)//'.struct',status='old')
   open(unit=unitkgen,file=clearspace(seedname)//'.outputkgen_orig',status='old')
   call countatoms(unitstruct,lattice)
   call readin_lattice(unitstruct,unitkgen,lattice)
   do jd1=1,3
     do jd2=1,3
        lattice%transform_matrix(jd1,jd2) = round_real(lattice%transform_matrix(jd1,jd2),4)
     enddo
   enddo


  !read-in klist
      open(klistunit,file=trim(seedname)//'.klist')
      nk=line_count(klistunit) - 2  !take care!!!!     
      !allocation of some arrays
      allocate(k(nk,3),kdiv(nk),wei(nk),k2(nk,3))
      do j=1,nk
        read(klistunit,"(A10,A80)")dum10,dumrest
        read(dumrest,*)k(j,:),kdiv(j)
!         if(j.eq.1) then 
          !           read(klistunit,1523) IK,(k(j,i),i=1,3),kdiv(j), &
!                         wei(j),emin,emax,nkp,kdim
          ! if (nk.ne.nkp) then
          !   write(*,*)"Error: number of k-points inconsistent, check .klist file"
  !           stop
  !        endif
!         else
!           read(klistunit,1520) IK, (k(j,i),i=1,3),kdiv(j),wei(j)
!         endif 
     enddo
     close(klistunit)
 1523 FORMAT(I10,4I10,3f5.1,4x,i6,10x,3i3,1x) 
 1520 FORMAT(I10,4I10,f5.1)
  
    k2 = 0d0
    do jk=1,nk
      
!       write(*,*)"tkpo1",k(jk,:)
      do j=1,3
          tmpvec(j) = real(k(jk,j))/real(kdiv(jk))
      enddo
!       write(*,*)"tkpo2",tmpvec
      if ((lattice%specifier.ne."R").and.(lattice%specifier.ne."H").and.(lattice%specifier.ne."C").and.(.not.lattconv)) then
            tmpvec = matmul(lattice%transform_matrix,tmpvec)
      endif
      k2(jk,:) = tmpvec
!       write(*,*)"tkpo3",k2(jk,:)
    enddo
    


  open(hamunit,file=trim(seedname)//'_hr.dat',status='old')
!   open(targethamunit,file=trim(seedname)//'.ham_fine',status='unknown')
  read(hamunit,*)dummy
  read(hamunit,*)nb
  read(hamunit,*)nr
  allocate(rweights(nr),Hr(nr,nb,nb),rvec(nr,3))
  allocate(Hk(nk,nb,nb),loc_ham(2*nb))
  rweights = 0
  do jr=1,nr/15
    read(hamunit,*)rweights((jr-1)*15+1:min(jr*15,nr))
!      write(*,"(16I4)")jr,rweights((jr-1)*15+1:min(jr*15,nr))
  enddo
  if (nr/15.ne.dble(nr)/15d0) then
      jr = nr/15+1
      read(hamunit,*)rweights((jr-1)*15+1:min(jr*15,nr))
  endif 
  Hr = dcmplx(0d0,0d0)
!    write(*,*)"nk,nr,nb",nk,nr,nb,Hr(1,1,1)
  do jr=1,nr
     do jb1=1,nb
        do jb2=1,nb
            read(hamunit,*)rvec(jr,:),di1,di2,realval,imagval
           Hr(jr,di1,di2) = dcmplx(realval,imagval)
!             write(*,*)rvec(jr,:),di1,di2,Hr(jr,di1,di2)
        enddo
     enddo
  enddo
  
  Hk = dcmplx(0d0)
  do jk=1,nk
      do jr=1,nr
             rdotk=2d0*pi*dot_product(k2(jk,:),rvec(jr,:))
             fac=exp(Xi*rdotk)/real(rweights(jr))
!              write(*,*)jk,jr
             do jb1=1,nb
                  do jb2=1,nb
                    Hk(jk,jb1,jb2)=Hk(jk,jb1,jb2)+fac*Hr(jr,jb1,jb2)
                 enddo
             enddo
      end do
  enddo
  
 !write-out
 open(unit=targethamunit,file=trim(seedname)//'.ham_fine',status='unknown')
 header=' # kpoints wann bands'
 write(targethamunit,'(I10,I6,I6,A21)')nk,nb,nb,header
! write(*,*)Hk(1,2,:)
! stop
 do jk=1,nk
     write(targethamunit,'(3F12.8)')k2(jk,:)
      do jb1=1,nb
          do jb2=1,nb
              loc_ham(2*jb2-1) = real(Hk(jk,jb1,jb2))
              loc_ham(2*jb2) = aimag(Hk(jk,jb1,jb2))
          enddo
          write(targethamunit,"(200F15.10)") loc_ham
!           write(*,"(30F15.10)")loc_ham
      enddo
!       stop
  enddo
  close(targethamunit)
  write(*,*)"Fourier-transformed Hamiltonian from _hr.dat to .ham_fine file according .klist file"

  END PROGRAM
