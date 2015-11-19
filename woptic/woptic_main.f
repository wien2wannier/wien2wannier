      program compute_optcond
      !main program to compute the optical conductivity via adaptive kmesh refinement

       use util 

      implicit none

        
      integer jj,iarg,matelmode,mommatunit,i1,i2,i3,i4,jk,n
      integer Nb,ii,jd,outputkgenunit,nt,ntfull
      integer :: klistunit,hamunit,outputunit,nemin,nemax,nkfull
      real*8 qxte,qyte,qzte,btransform(3,3)
      real*8 Pi,ediff,jointfilesum,sumrule,sumruledrude,sigmakfilesum,rvec(3)
      complex*16 Xi,cvec(3)
      integer, allocatable :: k(:,:),kperiodic(:,:),kmesh(:,:,:),kdiv(:),kfull(:,:)!,kw90(:,:)
      integer, allocatable :: tetra(:,:),tetrafull(:,:),map(:,:),rep(:,:),symopmap(:)
      integer, allocatable :: pmap(:),kidxt(:,:),kcandidate(:,:),kweights(:,:),patches(:,:,:)
      integer :: kdim(3),kvec(3)
      integer :: jd1,jd2,jd3,jd4,js,nop!,kvecw(3)
      integer nk,nkp,i,j,ik,woptinunit,rotvec(3)
      integer nbmax,found,internalhopping,jb1
      integer :: invidx,neminloc,nemaxloc,offset1,offset2
      real*8, allocatable :: qv(:,:),aa(:,:),wei(:),kcontrib(:),distmatrix(:,:,:)
      real*8, allocatable :: wtetra(:),symop(:,:,:),kdist(:,:,:)
      real*8, allocatable :: dxweight(:),dyweight(:),dzweight(:),bands(:,:)
      real*8 emin,emax,rbasis(3,3)
      real*8 krot(3)
      real*8 :: dkx,dky,dkz,dk(3),Vunit,dumsum(3)
      complex*16 :: Hloctet(10),derivative(3),locder(3),locval(3)
      integer :: foundk(10),nemin_wien2k,nemax_wien2k,nb_wien2k
      integer :: bandsunit,jb,nbloc,unitstruct,unitkgen,tmpunit
      complex*16, allocatable :: h(:,:),Hk(:,:,:),Hkdx(:,:,:),Hkfull(:,:,:)
      complex*16, allocatable :: Hkdy(:,:,:),Hkdz(:,:,:),mapcount(:)
      complex*16, allocatable :: Hkdxfull(:,:,:),Hkdyfull(:,:,:),Hkdzfull(:,:,:)
      complex*16, allocatable :: Hpolel(:,:,:)
      complex*16, allocatable :: polelx(:,:,:),ek(:,:,:),polely(:,:,:),polelz(:,:,:)
      parameter(Xi=(0.d0,1.d0))
      parameter(nbmax = 99)
      parameter (Pi = dacos(-1.d0))
      character*70 argdummy
      character*50 seedname
      character*120 dummy
      character*4 mode
      real*8 mxr,mxi,myr,myi,mzr,mzi,dw,sumr,drudesep,Fermienergy
      logical :: useinversion
      logical band
      character*10 dum10
  character*80 dumrest

      !logical :: WAN
      !character*3 :: wstr

      !definition for umat read-in
      integer, parameter :: dp2 = selected_real_kind(15,300)
      complex(kind=dp2), allocatable :: u_matrix(:,:,:)
       real*8 :: rotmatrix(3,3)
       integer :: jt!nemin_wann,nemax_wann,
       character*1 dum1
     type(structure) lattice

!     initialization
      hamunit = 11
      outputunit = 6
      woptinunit = 15
      klistunit = 8
      outputkgenunit = 9
      mommatunit = 12
      bandsunit = 16
      unitstruct = 17
      unitkgen = 18
      tmpunit = 19
      band = .false.
!       argument read-in
      iarg=iargc() 
      if(iarg.ge.1) then
         do j=1,iarg
            call getarg(j,argdummy)
            if (argdummy(1:1).eq.'-') then
                 if ((argdummy(2:2).eq.'h').or.(argdummy(2:2).eq.'H')) then
                  write(*,*)"compute optical conductivity for given k-mesh"
                  write(*,*)"input: case.woptin: input file"
                  write(*,*)"       case.klist: symmetrized k-points"
                  write(*,*)"       case.klist_full: unsymmetrized k-points"
                  write(*,*)"       case.kgen: symmetrized tetrahedra"
                  write(*,*)"       case.kgen_full: unsymmetrized tetrahedra"
                  write(*,*)"       case.kcontribw: contributions from previous iteration"
                  write(*,*)"                       '# 0 0 0' if first iteration or rerun"
                  write(*,*)"       case.map: internal mapping of klist_full to klist"
                  write(*,*)"       case.ham_fine: Hamiltonian H(k) in Wannier basis"
                  write(*,*)"                      refers to .klist_full for matelmode=1"
                  write(*,*)"                      refers to .klist for matelmode=2,3"
                  write(*,*)"       case.mommat: matrix elements for matelmode=2,3"
                  write(*,*)"output: case.woptoutw: standard output"
                  write(*,*)"        case.optcondw: optical conductivity"
                  write(*,*)"        case.kcontribw: contributions for this run"
                  write(*,*)"options: no flags, all options set in case.woptin"
                  stop
                elseif ((argdummy(2:5).eq.'band')) then
                    band = .true.
                else
                    write(*,*) 'Usage: woptic_main case [-band]'
                endif
            else!if (trim(argdummy).ne.'') then
               read(argdummy,*)seedname
            endif
         enddo
  
      else
         write(*,*) 'Usage: woptic_main case [-band]'
         stop
      endif
      !open output file
      open(outputunit,file=trim(seedname)//'.woptoutw',status='unknown')

       write(outputunit,*) "----- OPTICAL CONDUCTIVITY WITH WANNIER FUNCTIONS -----"

      !read-in input
      open(woptinunit,file=trim(seedname)//'.woptin',status='old')
      read(woptinunit,*)mode,matelmode
      read(woptinunit,*)dummy!emax,dw,deltino,beta
      read(woptinunit,*)nemin_wien2k,nemax_wien2k,nemin,nemax
      read(woptinunit,*)internalhopping
      read(woptinunit,*)drudesep
      close(woptinunit)
      nb = nemax - nemin + 1
      nb_wien2k = nemax_wien2k - nemin_wien2k + 1
      offset1= nemin - nemin_wien2k 
      offset2= nemax_wien2k - nemax 
      
     

      !read-in from seed.struct/seed.outputkgen
      open(unit=unitstruct,file=clearspace(seedname)//'.struct',status='old')
      open(unit=unitkgen,file=clearspace(seedname)//'.outputkgen_orig',status='old')
      call countatoms(unitstruct,lattice)
      call readin_lattice(unitstruct,unitkgen,lattice)
      close(unitstruct)
      close(unitkgen)

      !read-in klist
      if (band) then
          open(klistunit,file=trim(seedname)//'.klist_band')
      else
         open(klistunit,file=trim(seedname)//'.klist')
      endif
      nk=line_count(klistunit) - 2  !take care!!!!     
      !allocation of some arrays
      allocate(k(nk,3),kdiv(nk),wei(nk))
      allocate(Hk(nk,nb_wien2k,nb_wien2k),Hkdx(nk,nb_wien2k,nb_wien2k))
      Hk = dcmplx(0d0)
      Hkdx = dcmplx(0d0)
      allocate(polelx(nk,nb_wien2k,nb_wien2k),ek(nk,nb_wien2k,nb_wien2k),Hpolel(nk,nb_wien2k,nb_wien2k))
      ek = dcmplx(0d0)
      allocate(polely(nk,nb_wien2k,nb_wien2k),polelz(nk,nb_wien2k,nb_wien2k))
      allocate(Hkdy(nk,nb_wien2k,nb_wien2k),Hkdz(nk,nb_wien2k,nb_wien2k))
      Hkdy = dcmplx(0d0)
      Hkdz = dcmplx(0d0)
      allocate(aa(nb,2*nb),h(nb,nb),bands(nk,nb_wien2k))
      do j=1,nk
         read(klistunit,"(A10,A80)")dum10,dumrest
         read(dumrest,*)k(j,:),kdiv(j)
!         if(j.eq.1) then 
!           read(klistunit,1523) IK,(k(j,i),i=1,3),kdiv(j), &
!                         wei(j),emin,emax,nkp,kdim
!           if (nk.ne.nkp) then
!              write(outputunit,*)"Error: number of k-points inconsistent, check .klist file"
!              stop
!           endif
!         else
!           read(klistunit,1520) IK, (k(j,i),i=1,3),kdiv(j),wei(j)
!         endif 
     enddo
     close(klistunit)
 1523 FORMAT(I10,4I10,3f5.1,4x,i6,10x,3i3,1x) 
 1520 FORMAT(I10,4I10,f5.1)

      !read-in Fermi energy from scf2 file
      open(unit=tmpunit,file=clearspace(seedname)//'.fermienergy',status='old')
      read(tmpunit,"(A39,F9.5)")dummy,Fermienergy
      close(tmpunit)

     
      write(outputunit,*) "number of k points/ Wannier functions: ",nk,nb
      write(outputunit,*) "Fermi energy in wien2k is: ",Fermienergy,"Ry"
      write(outputunit,*) "number of elements: ",lattice%elementcount
      write(*,*)"lower and upper offset of energy windows: ",offset1,offset2 
!       if (nb.ne.nemax-nemin+1) then
!           write(outputunit,*)"Error: inconsistent number of Wannier orbitals, check .ham and .woptin"
!           stop
!       endif
      
     open(unit=44,file=trim(seedname)//".symop")
     read(44,*)nop
     allocate(symop(nop,3,3))
     do js=1,nop
        do jd=1,3
           read(44,"(3F8.5)")symop(js,jd,:)
!            write(*,*)js,lattice%sym(js,jd,:)
        enddo
        read(44,*)
     enddo
     close (44)

      !readin tetrahedra
      open(177,file=trim(seedname)//'.kgen')
      read(177,*)dum1,nt,kdim
      allocate(tetra(nt,10),wtetra(nt))
      allocate(kcandidate(kdim(1),4))
      do jt=1,nt
         read(177,*)tetra(jt,:),wtetra(jt)
      enddo
      close(177)

      !read-in information about unit cell and reciprocal basis
      open(outputkgenunit,file=trim(seedname)//'.outputkgen_orig',status='old')
      read(outputkgenunit,'(A120)')dummy
      do while (dummy(5:6).ne."G1") 
            read(outputkgenunit,'(A120)')dummy
      enddo
      !read over the reciprocal vectors without the 2pi
      read(outputkgenunit,'(A120)')dummy
      read(outputkgenunit,'(A120)')dummy
      read(outputkgenunit,'(A120)')dummy
      read(outputkgenunit,'(A120)')dummy
      write(outputunit,*)' reciprocal unit vectors in 2*pi*bohr^-1:'
      do j=1,3
         read(outputkgenunit,"(2x,3F10.6)")(rbasis(j,ii),ii=1,3)
         write(*,"(2x,3F10.6)"),rbasis(j,1:3)
      enddo

       

      !rotmatrix rotates a vector from the primitive basis to cartesian coordinates
      call inverse3x3(rbasis,rotmatrix)

      dkx = 0d0;
      dky = 0d0;
      dkz = 0d0;
      do j=1,3
         dkx = dkx + rbasis(j,1)**2;
         dky = dky + rbasis(j,2)**2;
         dkz = dkz + rbasis(j,3)**2;
      enddo
      dkx = 1d0
      dky = 1d0
      dkz = 1d0
      !write(*,*)dsqrt(dkx),dsqrt(dky),dsqrt(dkz),kdim
      !dk(1) = dsqrt(dkx)/dble(kdim(1))
      !dk(2) = dsqrt(dky)/dble(kdim(2))
      !dk(3) = dsqrt(dkz)/dble(kdim(3))
      dk(1) = 1d0/dble(kdim(1))
      dk(2) = 1d0/dble(kdim(2))
      dk(3) = 1d0/dble(kdim(3))
      
      write(outputunit,*)' reciprocal increments are [dkx dky dkz]:',dk(1),dk(2),dk(3)
      !compute unit cell volume in bohr^3
      Vunit = rbasis(1,1)*(rbasis(2,2)*rbasis(3,3)-rbasis(3,2)*rbasis(2,3)) &
              + rbasis(1,2)*(rbasis(2,3)*rbasis(3,1)-rbasis(3,3)*rbasis(2,1)) &
              + rbasis(1,3)*(rbasis(2,1)*rbasis(3,2)-rbasis(3,1)*rbasis(2,2)) 
      Vunit = 1d0/Vunit*(pi*2d0)**3
      write(outputunit,*)' unit cell volume in bohr^3:',Vunit
      !test values for SrVO3
      !dkx = 0.865298/real(kdim(1));
      !dky = 0.865298/real(kdim(2));
      !dkz = 0.865298/real(kdim(3));
      close(outputkgenunit)

      if (((offset1.ne.0).or.(offset2.ne.0)).or.(matelmode.eq.4)) then
!       readin bands
       write(outputunit,*)'>>> bands read in from .energy'
                open(bandsunit,file=trim(seedname)//'.energy',status='old')
	        do j=1,lattice%elementcount
			read(bandsunit,*)
			read(bandsunit,*)
		enddo
		do jk=1,nk
! 			read(bandsunit,"(3e19.12,i10,2i6,f5.1)")mxr,mxi,myr,ii,jj,nbloc
                        read(bandsunit,"(73X,i6)")nbloc
			write(*,*)nbloc
			do jb=1,nbloc
				read(bandsunit,*)ii,mxr
				flush(6)
				if ((jb.lt.nemin_wien2k).or.(jb.gt.nemax_wien2k)) then
					cycle
				else
					bands(jk,jb-nemin_wien2k+1) = (mxr-Fermienergy)*13.60514d0
					write(*,*)jk,jb,bands(jk,jb-nemin_wien2k+1)
				endif
			enddo
		enddo
		close(bandsunit)
      endif
        
     

      Hkdx(:,:,:) = dcmplx(0d0,0d0)
      Hkdy(:,:,:) = dcmplx(0d0,0d0)
      Hkdz(:,:,:) = dcmplx(0d0,0d0)

!       !some comparison readin
!       open(201,file=trim(seedname)//'.joint')
!       n=line_count(201) !take care!!!!     
!       read(201,*)
!       read(201,*)
!       read(201,*)
!       sumr = 0d0
!       do j=1,n-3
!          read(201,*)mxr,mxi,myr
!          if (j.eq.2)then
!            dw = mxr
!          endif
!          sumr = sumr + mxi
!       enddo
!       jointfilesum = sumr*dw
!       write(*,*)"joint file sum:",jointfilesum
!       close(201)
!       open(201,file=trim(seedname)//'.sigmak')
!       n=line_count(201) !take care!!!!     
!       do j=1,9
!             read(201,*)
!       enddo
!       sumr = 0d0
!       dw = 0d0
!       do j=1,n-9
!          read(201,*)mxr,mxi,myr
! !          if (j.eq.1).and.(mxr.eq.0d0) then
! !            dw = mxr
!          if (j.le.2) then
!            dw = mxr - dw
!          endif
!          if ((j.le.2).or.(dw*dble(j).le.drudesep)) then
!             sumr = sumr + mxi
!          endif
!       enddo
!       sigmakfilesum = sumr*dw
!       write(*,*)"sigmak file Drude sum,dw:",sigmakfilesum,dw
!       close(201)
      if (matelmode.ge.2) then
         write(outputunit,*)'>>> use polarization matrix elements'

!           nemin = 12
          open(mommatunit,file=trim(seedname)//'.mommat',status='old')
          read(mommatunit,"(A120)")dummy

          do j=1,nk
             read(mommatunit,"(A120)")dummy
             read(mommatunit,"(A120)")dummy
             read(dummy(28:38),"(2I5)")neminloc,nemaxloc
             read(mommatunit,"(A120)")dummy
             do ii=neminloc,nemaxloc
               do jj=ii,nemaxloc
                  if (jj.ge.ii) then
                       read(mommatunit,"(A120)")dummy
                       read(dummy,9040)i1,i2,mxr,mxi,myr,myi,mzr,mzi,ediff
                       if ((ii.lt.nemin_wien2k).or.(ii.gt.nemax_wien2k).or.  &
                           (jj.lt.nemin_wien2k).or.(jj.gt.nemax_wien2k)) then
                            cycle
                       else
                           Hkdx(j,ii-nemin_wien2k+1,jj-nemin_wien2k+1) = dcmplx(mxr,mxi)
                           Hkdy(j,ii-nemin_wien2k+1,jj-nemin_wien2k+1) = dcmplx(myr,myi)
                           Hkdz(j,ii-nemin_wien2k+1,jj-nemin_wien2k+1) = dcmplx(mzr,mzi)
                       endif
                       if (jj.ne.ii) then
                         Hkdx(j,jj-nemin_wien2k+1,ii-nemin_wien2k+1) = &
                                        conjg(Hkdx(j,ii-nemin_wien2k+1,jj-nemin_wien2k+1))
                         Hkdy(j,jj-nemin_wien2k+1,ii-nemin_wien2k+1) = &
                                         conjg(Hkdy(j,ii-nemin_wien2k+1,jj-nemin_wien2k+1))
                         Hkdz(j,jj-nemin_wien2k+1,ii-nemin_wien2k+1) = &
                                         conjg(Hkdz(j,ii-nemin_wien2k+1,jj-nemin_wien2k+1))
                       endif
                  endif
               enddo
             enddo
          enddo
     endif 

     !obtain groupvelocities according to option in woptin file
     if (matelmode.eq.1) then !Peierls

      !we need the full kmesh only in the Peiersl case
      open(klistunit,file=trim(seedname)//'.klist_full')
      nkfull = line_count(klistunit)-2 !take care!!!! 
      !write(*,*)"loner",line_count(klistunit),klistunit
      allocate(kfull(nkfull,3),kmesh(3,nkfull,2),kdist(3,nkfull,2))
      allocate(Hkdxfull(nkfull,nb_wien2k,nb_wien2k),Hkdyfull(nkfull,nb_wien2k,nb_wien2k),Hkdzfull(nkfull,nb_wien2k,nb_wien2k))
      allocate(map(nkfull,2),mapcount(nk),dxweight(nkfull),dyweight(nkfull),dzweight(nkfull),Hkfull(nkfull,nb_wien2k,nb_wien2k))
      Hkfull = dcmplx(0d0)
      allocate(rep(nk,2),symopmap(nkfull),pmap(nkfull),kidxt(nkfull,3),kweights(nkfull,3))
     
      flush(6)
      do j=1,nkfull
        if(j.eq.1) then 
          read(klistunit,1523) IK,(kfull(j,i),i=1,3),ii, &
                        mxr,emin,emax,nkp,kvec
          if (nkfull.ne.nkp) then
             write(outputunit,*)"Error: number of k-points inconsistent, check .klist_full file"
             stop
          endif
        else
          read(klistunit,1520) IK, (kfull(j,i),i=1,3),ii,mxr
        endif 
      enddo
      close(klistunit)
      do j=1,nkfull
         !write(*,*)"kfull before map",j,lattice%transform_matrix(1,:),kfull(j,:)
         do jd1=1,3
            rvec(jd1) = real(kfull(j,jd1))/real(kdim(jd1))
         enddo
         if ((lattice%specifier.ne."R").and.(lattice%specifier.ne."H")) then
            rvec = matmul(lattice%transform_matrix,rvec)
         endif
         do jd1=1,3
            kfull(j,jd1) = int(rvec(jd1)*kdim(jd1))
         enddo
         !write(*,*)"kfull after map",j,kfull(j,:)
      enddo
      !write(*,*)"kfull2",kfull(1:3,:)
      call inverse3x3(lattice%transform_matrix,btransform)
      !write(*,*)"bvec1",btransform(1,:)
      !write(*,*)"bvec2",btransform(2,:)
      !write(*,*)"bvec3",btransform(3,:)
 

      open(177,file=trim(seedname)//'.kgen_full')
      read(177,*)dum1,ntfull,kdim
      allocate(tetrafull(ntfull,10))
      do jt=1,ntfull
         read(177,*)tetrafull(jt,:),mxr
!               write(*,*)"tf",tetrafull(jt,:),mxr
      enddo
      close(177)
      open(177,file=trim(seedname)//'.map')
      do jk=1,nkfull
         read(177,*)map(jk,:)
      enddo
      close(177)

      !in case of Peierls the full Hamilonian without symmetries is read in
      open(hamunit,file=trim(seedname)//'.ham_fine',status='old')
      read(hamunit,*)  nkp,nb
      if (nkfull.ne.nkp) then
          write(outputunit,*)"Error: number of k-points inconsistent, check .klist_full and .ham_fine file"
          stop
      endif
      allocate(qv(nkfull,3),patches(nkfull,3072,3))
!     read-in Hamiltonian
      do i=1,nkfull
        read(hamunit,*) qxte,qyte,qzte
        !write(*,*)qxte,qyte,qzte
        qv(i,1)=qxte
        qv(i,2)=qyte
        qv(i,3)=qzte
        do j=1, Nb
           read(hamunit,*) (aa(j,ii),ii=1,2*Nb)
        enddo
        do j=1, Nb
            do ii=1,Nb
               Hkfull(i,offset1+j,offset1+ii)=aa(j,2*ii-1)+Xi*aa(j,2*ii)
!                write(*,*)i,j,ii,Hkfull(i,offset1+j,offset1+ii)
            enddo
        enddo
     enddo
     close(hamunit)

     !now fill up the Hamiltonian with band energies if required
     do i=1,nkfull
            do j=1, offset1
!                write(*,*)"written-",j
               if (bands(map(i,1),j).eq.0d0) then
                    write(outputunit,*)"Warning: lower index of bands is lower than overall wien2k index"
               endif
               Hkfull(i,j,j) = bands(map(i,1),j)
            enddo
            do j=1, offset2
!                 write(*,*)"written+",nemax-nemin_wien2k+j+1
               if (bands(map(i,1),nemax-nemin_wien2k+j+1).eq.0d0) then
                    write(outputunit,*)"Warning: upper index of bands is larger than overall wien2k index"
               endif
               Hkfull(i,nemax-nemin_wien2k+j+1,nemax-nemin_wien2k+j+1) = bands(map(i,1),nemax-nemin_wien2k+j+1)
           enddo
        enddo    
     write(outputunit,*)'>>> finished reading Hamiltonian'
     write(outputunit,*)'>>> compute group velocity via derivative of Hamiltonian ...'
     !compute the derivative on the full kmesh
      Hkdxfull = cmplx(0d0)
      Hkdyfull = cmplx(0d0)
      Hkdzfull = cmplx(0d0)
      patches = 0
      do jt=1,ntfull
          do jd1=1,10
              do jd2=1,nb_wien2k
                  do jd3=1,nb_wien2k
                       Hloctet = Hkfull(tetrafull(jt,:),jd2,jd3)
                       call compute_tetrahedralderivative(nkfull,tetrafull(jt,:),&
                                            kfull,dk,Hloctet,jd1,kdim,derivative,0,rbasis)
                       cvec = derivative!matmul(dcmplx(rotbasis),derivative)
                       Hkdxfull(tetrafull(jt,jd1),jd2,jd3) = Hkdxfull(tetrafull(jt,jd1),jd2,jd3) + cvec(1)
                       Hkdyfull(tetrafull(jt,jd1),jd2,jd3) = Hkdyfull(tetrafull(jt,jd1),jd2,jd3) + cvec(2)
                       Hkdzfull(tetrafull(jt,jd1),jd2,jd3) = Hkdzfull(tetrafull(jt,jd1),jd2,jd3) + cvec(3)
                  enddo
              enddo
              patches(tetrafull(jt,jd1),1,1) = patches(tetrafull(jt,jd1),1,1) + 1
          enddo
      enddo
      
      !normalize
      do jk=1,nkfull
         do jd1=1,nb_wien2k
             do jd2=1,nb_wien2k
                  Hkdxfull(jk,jd1,jd2) = Hkdxfull(jk,jd1,jd2)/cmplx(patches(jk,1,1))
                  Hkdyfull(jk,jd1,jd2) = Hkdyfull(jk,jd1,jd2)/cmplx(patches(jk,1,1))
                  Hkdzfull(jk,jd1,jd2) = Hkdzfull(jk,jd1,jd2)/cmplx(patches(jk,1,1))
             enddo
         enddo
      enddo

      !now map to symmetrized mesh
      mapcount = 0
      do jk=1,nkfull
         if (mapcount(map(jk,1)).eq.0) then
		do jd1=1,nb_wien2k
		do jd2=1,nb_wien2k
			locder(1) = Hkdxfull(jk,jd1,jd2)
			locder(2) = Hkdyfull(jk,jd1,jd2)
			locder(3) = Hkdzfull(jk,jd1,jd2)
			do jd3=1,3
                                Hkdx(map(jk,1),jd1,jd2) = Hkdx(map(jk,1),jd1,jd2) + symop(map(jk,2),1,jd3)*locder(jd3)
				Hkdy(map(jk,1),jd1,jd2) = Hkdy(map(jk,1),jd1,jd2) + symop(map(jk,2),2,jd3)*locder(jd3)
				Hkdz(map(jk,1),jd1,jd2) = Hkdz(map(jk,1),jd1,jd2) + symop(map(jk,2),3,jd3)*locder(jd3)
			enddo
			Hk(map(jk,1),jd1,jd2) = Hkfull(jk,jd1,jd2)
		enddo
		enddo
		
		mapcount(map(jk,1)) = mapcount(map(jk,1)) + 1 
          endif
      enddo

      !normalize
      do jk=1,nk
         !write(*,*)jk,k(jk,:)
         do jd1=1,nb_wien2k
             do jd2=1,nb_wien2k
                  Hkdx(jk,jd1,jd2) = Hkdx(jk,jd1,jd2)/mapcount(jk)
                  Hkdy(jk,jd1,jd2) = Hkdy(jk,jd1,jd2)/mapcount(jk)  
                  Hkdz(jk,jd1,jd2) = Hkdz(jk,jd1,jd2)/mapcount(jk) 
             enddo
         enddo
      enddo
   
    !internal hopping
     if (internalhopping.eq.1) then
         write(outputunit,*)'>>> add internal unitcell hopping...'
         !read-in distance matrix
         allocate(distmatrix(3,nb,nb))
         open(unit=177,file=trim(seedname)//'.distmatrix',status='old')
         do jd1=1,3
            do jb1=1,nb
                read(177,*)distmatrix(jd1,jb1,:)
            enddo
         enddo
         close(177)
         do jk=1,nk
               do jd1=1,nb
                  do jd2=1,nb
                      Hkdx(jk,offset1+jd1,offset1+jd2) = Hkdx(jk,offset1+jd1,offset1+jd2)&
                  - cmplx(0d0,1d0)*distmatrix(1,jd1,jd2)*Hk(jk,jd1,jd2)
                      Hkdy(jk,offset1+jd1,offset1+jd2) = Hkdy(jk,offset1+jd1,offset1+jd2)&
                  - cmplx(0d0,1d0)*distmatrix(2,jd1,jd2)*Hk(jk,jd1,jd2)
                      Hkdz(jk,offset1+jd1,offset1+jd2) = Hkdz(jk,offset1+jd1,offset1+jd2)&
                  - cmplx(0d0,1d0)*distmatrix(3,jd1,jd2)*Hk(jk,jd1,jd2)
                  enddo
               enddo
         enddo
     endif
   
    !do jk=1,nk
         !write(*,*)jk,k(jk,:)
         !do jd1=1,nb_wien2k
         !    do jd2=1,nb_wien2k
         !          write(195,*)jk,jd1,jd2,Hk(jk,jd1,jd2),Hkdx(jk,jd1,jd2),Hkdy(jk,jd1,jd2),Hkdz(jk,jd1,jd2) 
         !    enddo
         !enddo
     ! enddo

         

       elseif (matelmode.ge.2) then!read dipole matrix elements from wien2k

           if ((matelmode.eq.2).or.(matelmode.eq.3)) then
              allocate(u_matrix(nb_wien2k,nb_wien2k,nk))
!            do i=1,nb_wien2k
!              do j=1,nb_wien2k
! !                  write(*,*)"Hk0a",i,j,Hk(2,i,j)
!              enddo
!          enddo
               
             !in case of dipole matrix elements the symmetrized Hamiltonian is readin
		open(hamunit,file=trim(seedname)//'.ham_fine',status='old')
		read(hamunit,*)  nkp,nb
                allocate(qv(nk,3))
		!     read-in Hamiltonian
		do i=1,nk
			read(hamunit,*) qxte,qyte,qzte
			qv(i,1)=qxte
			qv(i,2)=qyte
			qv(i,3)=qzte
			do j=1, Nb
			read(hamunit,*) (aa(j,ii),ii=1,2*Nb)
			enddo
			do j=1, Nb
			do ii=1,Nb
			Hk(i,offset1+j,offset1+ii)=aa(j,2*ii-1)+Xi*aa(j,2*ii)
			enddo
			enddo
		enddo
		close(hamunit)
         !write(outputunit,*)"offset1,offset2",offset1,offset2,nemin_wien2k,nemax_wien2k,nemin,nemax
!           do i=1,nb_wien2k
!              do j=1,nb_wien2k
!                  write(*,*)"Hk0",i,j,Hk(2,i,j)
!              enddo
!          enddo
         do i=1,nk
            do j=1, offset1
!                write(*,*)"written-",j
               if (bands(i,j).eq.0d0) then
                    write(outputunit,*)"Warning: lower index of bands is lower than overall wien2k index"
               endif
               Hk(i,j,j) = bands(i,j)
            enddo
            do j=1, offset2
!                 write(*,*)"written+",nemax-nemin_wien2k+j+1
               if (bands(i,nemax-nemin_wien2k+j+1).eq.0d0) then
                    write(outputunit,*)"Warning: upper index of bands is larger than overall wien2k index"
               endif
               Hk(i,nemax-nemin_wien2k+j+1,nemax-nemin_wien2k+j+1) = bands(i,nemax-nemin_wien2k+j+1)
           enddo
        enddo            
      !do i=1,nb_wien2k
      !    do j=1,nb_wien2k
      !       write(*,*)"Hk1",i,j,Hk(1,i,j)
      !    enddo
      !enddo
         
        do i=1,nk
            call get_umatrix(Hk(i,1:nb_wien2k,1:nb_wien2k),nb_wien2k,u_matrix(1:nb_wien2k,1:nb_wien2k,i))  
        enddo


           !rotate Hamiltonian back to diagonal basis
           do jk=1,nk
             do i1=1,nb_wien2k
               do i2=1,nb_wien2k
                  do i3=1,nb_wien2k
                     do i4=1,nb_wien2k
                         ek(jk,i1,i2) = ek(jk,i1,i2) + &
               conjg(u_matrix(i3,i1,jk))*Hk(jk,i3,i4)*u_matrix(i4,i2,jk)
                     enddo
                  enddo
               enddo
             enddo
          enddo
!           stop
       
          !read-in matrix elements from mommat file
!           if ((matelmode.eq.2) .or. (matelmode.eq.3)) then
          endif    
       
            
          if (matelmode.eq.2) then!use rotated matrix elements
             write(outputunit,*)'  use rotated matrix elements'
             write(outputunit,*)'  with wannier90 Hamiltonian'

             flush(outputunit)
             do j=1,nk
               do ii=1,nb_wien2k
                 do jj=1,nb_wien2k          
                    do i1=1,nb_wien2k
                      do i2=1,nb_wien2k
                         polelx(j,ii,jj) = polelx(j,ii,jj) + &
                                        u_matrix(ii,i1,j)*Hkdx(j,i1,i2)*conjg(u_matrix(jj,i2,j))
                         polely(j,ii,jj) = polely(j,ii,jj) + &
                                        u_matrix(ii,i1,j)*Hkdy(j,i1,i2)*conjg(u_matrix(jj,i2,j))
                         polelz(j,ii,jj) = polelz(j,ii,jj) + &
                                        u_matrix(ii,i1,j)*Hkdz(j,i1,i2)*conjg(u_matrix(jj,i2,j))
                      enddo
                    enddo
                 enddo
              enddo
              do ii=1,nb_wien2k
                 do jj=1,nb_wien2k   
                     Hkdx(j,ii,jj) = polelx(j,ii,jj)
                     Hkdy(j,ii,jj) = polely(j,ii,jj)
                     Hkdz(j,ii,jj) = polelz(j,ii,jj)
                 enddo
              enddo
            enddo
            do j=1,nk
                do ii=1,nb_wien2k
!                   do jj=1,nb_wien2k
!                       Hk(j,ii,jj) = ek(j,ii,jj)
                     
!                        Hk(j,ii,ii) = bands(j,ii) 
!                   enddo
                  write(*,"(40F12.5)") Hk(j,ii,:)
                enddo
             enddo
             write(*,*)"interlagos"
             do j=1,nk
                do ii=1,nb_wien2k
!                   do jj=1,nb_wien2k
!                       Hk(j,ii,jj) = ek(j,ii,jj)
                     
!                        Hk(j,ii,ii) = bands(j,ii) 
!                   enddo
                  write(*,"(40F12.5)") u_matrix(ii,:,j)
                enddo
             enddo
!              stop

          elseif (matelmode.eq.3) then!use diagonal Hamiltonian
            write(outputunit,*)'  use original wien2k matrix elements'
            write(outputunit,*)'  with diagonal Hamiltonian'
!            Hk = dcmplx(0d0)
             do j=1,nk
                do ii=1,nb_wien2k
                  do jj=1,nb_wien2k
                      Hk(j,ii,jj) = ek(j,ii,jj)
                     
!                        Hk(j,ii,ii) = bands(j,ii) 
                  enddo
                  write(*,"(40F12.5)") Hk(j,ii,:)
                enddo
             enddo
!             write(*,*)"running debugging function"
!              call optic_direct_tetra(seedname,nk,nb_wien2k,Hkdx,Hkdy,Hkdz&
!                                ,ek,wei,nt,tetra,wtetra)
!              stop
          elseif (matelmode.eq.4) then!wien2k only mode
             write(outputunit,*)'  use original wien2k matrix elements'
             write(outputunit,*)'  with diagonal Hamiltonian from wien2k'
             do i=1,nk
                do j=1, nb_wien2k
!                    if (bands(i,j).eq.0d0) then
!                        write(outputunit,*)"Warning: lower index of bands is lower than overall wien2k index"
!                   endif
                  Hk(i,j,j) = bands(i,j)
               enddo
!             do j=1, offset2
!                 write(*,*)"written+",nemax-nemin_wien2k+j+1
!                if (bands(i,nemax-nemin_wien2k+j+1).eq.0d0) then
!                     write(outputunit,*)"Warning: upper index of bands is larger than overall wien2k index"
!                endif
!                Hk(i,nemax-nemin_wien2k+j+1,nemax-nemin_wien2k+j+1) = bands(i,nemax-nemin_wien2k+j+1)
            enddo
!         enddo 

          else
             write(outputunit,*)"Error: did not recognize matrix element mode."
             write(outputunit,*)"check input file" 
             stop
          endif
           
          close(mommatunit)
       
          
       else
         write(outputunit,*)"Error: did not recognize matrix element mode."
         write(outputunit,*)"check input file" 
         stop
       endif
 9040 FORMAT(3X,2I4,6E13.6,F13.8)
 !9050 FORMAT(3X,2I4,6E19.11)
       
        write(outputunit,*)'>>> compute optical conductivity...'
        flush(6)
!        call debugfunct3(seedname,nk,k,nt,tetra,wtetra,lattice,kdim)    
        call optcalc(lattice,Vunit,Nb_wien2k,nk,nt,tetra,wtetra,Hk,Hkdx,Hkdy,Hkdz&
                         ,sumruledrude,sumrule,seedname,woptinunit,outputunit,band)
       write(*,*)'END, all clear'
       close(outputunit)
       
          
       CONTAINS
       FUNCTION kdiff(kp1,kp2,ndim) result(outm)
    
       implicit none

       integer kp1(3),kp2(3),ndim(3),outm(3),jd,kp1tmp(3),kp2tmp(3)
         
       outm = 0
       kp1tmp = kp1
       kp2tmp = kp2
       do jd=1,3
!            write(*,*)"yuppy",kp1(:),kp2(:),(kp1(jd).eq.0),(kp2(jd).gt.ndim(jd)/2)
           if ((kp1(jd).eq.0).and.(kp2(jd).gt.ndim(jd)/2)) then
!                write(*,*)"warning",kp1(:),kp2(:)
                kp1tmp(jd) = ndim(jd)
           endif
           if ((kp2(jd).eq.0).and.(kp1(jd).gt.ndim(jd)/2)) then
!                write(*,*)"warning",kp1(:),kp2(:)
                kp2tmp(jd) = ndim(jd)
           endif
!           if (kp2(jd).gt.ndim(jd)/2) then
!                kp2(jd) = ndim(jd)-kp2(jd)
!           endif
          outm(jd) = kp2tmp(jd)-kp1tmp(jd)
       enddo
       END FUNCTION

       end 
 
!     *****************************************************************
      FUNCTION find_kpoint(nk,k,kvec,kdim) result(outp)

      implicit none

      integer outp,nk,k(nk,3),kvec(3),j,kdim(3)

      outp = 0
      
      do j=1,nk
         if (k(j,1).gt.kdim(1)/2) then
              k(j,1) = kdim(1) - k(j,1) 
         endif
         if (k(j,2).gt.kdim(2)/2) then
            k(j,2) = kdim(2) - k(j,2)
         endif
         if (k(j,3).gt.kdim(3)/2) then
           k(j,3) = kdim(3) - k(j,3) 
          endif
           if (((k(j,1).eq.kvec(1)).or.(k(j,1)-kdim(1).eq.kvec(1)))&
          .and.((k(j,2).eq.kvec(2)).or.(k(j,2)-kdim(2).eq.kvec(2)))&
          .and.((k(j,3).eq.kvec(3)).or.(k(j,3)-kdim(3).eq.kvec(3)))) then
              outp = j
              exit
           endif
      enddo


      END FUNCTION

!     *****************************************************************
      SUBROUTINE optcalc(lattice,Vunit,Nb,nk,nt,tetra,wtetra,Hk,Hkdxfin,Hkdyfin,Hkdzfin,&
                     sumruledrude1,sumrule1,seedname,inputunit,outputunit,band)
      use util

      implicit none
      integer :: Nb,Nk,Nkold!,j1,j2,i1,i2
      integer :: i,il,jj,k,jd
      integer :: jb1,jb2,jk,jw,jd1,jd2,jd3,jd4,nsig,jsig,nop,reftet(8,4),jda
      integer :: maxw,inputunit,outputunit,idx1,idx2,jwext,jwint
      logical :: opt,sfwarn,nonint,compdos,diag,joint,resolvecont
      integer :: nt,tetra(nt,10),counter,js,jt,matelmode,intidx(40),ninter,SEnfreq
      integer, allocatable :: sigidx(:,:)
      real*8 :: beta,sumruledrude1,sumrule1
      real*8 :: pi,wtetra(nt)
      real*8 :: f_j,f_kj,emax, dw, deltino 
      real*8,allocatable ::  w(:),DOS(:,:,:),DOStot(:,:),totDOS(:),optcond(:,:,:,:),tmpoptcond(:,:,:)
      real*8,allocatable ::  totoptcond(:,:,:),symop(:,:,:),tmp(:),tmpDOS(:)
      real*8,allocatable ::  tmpr(:),optcond_resolved(:,:,:,:,:),totoptcond_resolved(:,:,:,:)
      real*8 :: sumtetraweights,df,totdccond(3,3),totdccond_resolved(3,3,nb),Vunit
      real*8, allocatable :: sumruledrude(:),sumrule(:),wSE(:),tmpSE(:)
      complex*16 :: Hk(Nk,Nb,Nb),Hkdxfin(Nk,Nb,Nb)!,umatrix(nb,nb,nk)
      complex*16 :: Hkdyfin(Nk,Nb,Nb),Hkdzfin(Nk,Nb,Nb)
      complex*16 :: Xi, eye(nb,nb)
      complex*16 :: LHS(Nb,Nb),RHS1(Nb,Nb),LHSc(Nb,Nb),RHS1c(Nb,Nb)
      complex*16 :: RHS2(Nb,Nb),RHS2c(Nb,Nb),matel1(nb,nb),matel2(nb,nb)
      complex*16, allocatable :: Gtvk1(:,:,:),Gtvkc1(:,:,:)
      complex*16, allocatable :: Gtvk2(:,:,:),Gtvkc2(:,:,:)
      complex*16, allocatable ::  SE(:,:),SEin(:,:)
      real*8 :: dccond(nk,3,3),convfac,echarge,hbar,K1(nk,3,3),totK1(3,3),thermopower(3,3),dwSE
      real*8 :: dccond_resolved(nk,3,3,nb),K1_resolved(nk,3,3,nb),totK1_resolved(3,3,nb),thermopower_resolved(3,3,nb)
      character*4 :: mode
      character*50 :: seedname
      character*1 dum1
      real*8 :: drudesep,factor1,minimum
      parameter (Xi=(0d0,1d0))
      parameter (pi=dacos(-1d0))
      type(structure) lattice
      logical :: applyscissors,foundw,band
      real*8 scissorsval
      integer :: scissorsnorb
      integer, allocatable :: scissorsorbidx(:)
      real*8 :: chempot

      !debug
!       Hkdxfin = dcmplx(1d0)
!       Hkdyfin = dcmplx(1d0)
!       Hkdzfin = dcmplx(1d0)

      resolvecont = .false.
!      write(*,*)"strange",Hkdxfin(1,1,10),Hkdyfin(1,1,10),Hkdzfin(1,1,10)
!      Hkdzfin = Hkdxfin !!!Attenzione!!!!
      !readin input
      open(inputunit,file=trim(seedname)//'.woptin',status='old')
      read(inputunit,*)mode,matelmode
!       write(*,*)mode,matelmode
      read(inputunit,*)emax,dw,deltino,beta
!       write(*,*)emax,dw,deltino,beta
      read(inputunit,*)dum1
!       write(*,*)dum1
      read(inputunit,*)dum1
!       write(*,*)dum1
      read(inputunit,*)drudesep
!       write(*,*)drudesep
      read(inputunit,*)nonint,resolvecont
      read(inputunit,*)chempot
      write(outputunit,*)"chemical potential is: ",chempot
      read(inputunit,*)applyscissors,scissorsval,scissorsnorb
      if (applyscissors) then
         allocate(scissorsorbidx(scissorsnorb))
         read(inputunit,*)scissorsorbidx
         write(outputunit,*)'energy shift of',scissorsval
         write(outputunit,*)'will be applied to the orbitals',scissorsorbidx
      endif
!       if (.not.noninit) then
!          read(inputunit,*)intidx
!       endif
      close(inputunit)

      !Vunit should be bohr^3
      echarge = 1.602176487d-19;![As]
      hbar = 1.05457148d-34;!Js
      factor1 = 1.413351709413265d+08!hbar^3/emass^2*(m/ang)^5
      write(*,*)"mode: ",matelmode
      write(*,*)"number of orbitals contributing to optics: ",nb
      if (matelmode.eq.1) then!Peierls
         !Expl:    2*pi*elem.charge^2/hbar/unitcellvol[bohr^3] * ang/bohr*meter/ang*cm/meter (*sigma[bohr^2] from Ham.)
         convfac = 2d0*pi*echarge**2/hbar/Vunit*(1d0/0.529177d0)*10E10/10E2!conv. to(Siemens/cm)
      else!wien2k matrix elements
         !Expl: echarge*hbar^3[Js]*(J/eV)^2/emass^2/bohr^5*(bohr^5/0.52^5 ang^5)*(10^50 ang^5/m^5)*(m/10^2 cm)
         !strange wien2k factor coming from joint.f: 64/pi???
         !spin is already included in matrix elements
         convfac = 64d0*factor1/Vunit/(0.529177**5)/10E2!conv. to(Siemens/cm)
      endif
      write(*,*)"convfac: ",convfac

      !initialization
      compdos = .false.
      opt=.false.
!       nonint = .true.
      joint = .false.
      if (band) then
        write(outputunit,*)'BAND mode: k-wise contributions'
        opt = .true.
        compdos = .true.
      elseif (trim(mode).eq."DOS") then
        compdos = .true.
        write(outputunit,*)'DOS mode: compute just the DOS'
      elseif (trim(mode).eq."OPT") then
        write(outputunit,*)'OPT mode: compute non-interacting opt. conductivity and DOS'
        opt = .true.
        compdos = .true.
      elseif (trim(mode).eq."IOPT") then
        write(outputunit,*)'IOPT mode: compute interacting opt. conductivity with'
        write(outputunit,*)'self energies from Sf_<bandidx>.dat and DOS'
        opt = .true.
        compdos = .true.
        nonint = .false.
      elseif (trim(mode).eq."JOIN") then
        write(outputunit,*)'JOIN modus: compute joint density of states'
        opt = .true.
        joint = .true.
        compdos = .true.
        Hkdxfin = dcmplx(1d0)
        Hkdyfin = dcmplx(1d0)
        Hkdzfin = dcmplx(1d0)
        convfac = 1d0
      endif

      maxw = 2*int(emax/dw)
      nsig = 6
      
      allocate(Gtvk1(-maxw:maxw,Nb,Nb),Gtvkc1(-maxw:maxw,Nb,Nb))
      allocate(Gtvk2(-maxw:maxw,Nb,Nb),Gtvkc2(-maxw:maxw,Nb,Nb))
      allocate(SE(-maxw:maxw,Nb))
      allocate(DOS(nk,-maxw:maxw,Nb),DOStot(-maxw:maxw,Nb),totDOS(-maxw:maxw),tmpDOS(nb))
      allocate(w(-maxw:maxw))
      allocate(optcond(nk,1:maxw-1,3,3),tmpoptcond(1:maxw-1,3,3),totoptcond(1:maxw-1,3,3))
      allocate(sigidx(nsig,2),tmp(nsig))
      sigidx(1,:) = (/1,1/)
      sigidx(2,:) = (/1,2/)
      sigidx(3,:) = (/1,3/)
      sigidx(4,:) = (/2,2/)
      sigidx(5,:) = (/2,3/)
      sigidx(6,:) = (/3,3/)
      allocate(sumruledrude(nsig))
      allocate(sumrule(nsig))

      dccond = 0d0
      K1 = 0d0
      optcond = 0d0
      sfwarn=.false.
      if (compdos) then
        DOS = 0d0
        DOStot = 0d0
        totDOS = 0d0
      endif
      write(outputunit,*)'parameters: [emax dw deltino beta]'
      write(outputunit,"(2x,4F13.6)")emax,dw,deltino,beta

      open (50,file=trim(seedname)//'.optcondw',status='unknown')
      open (51,file=trim(seedname)//'.wdosw',status='unknown')
      open (52,file=trim(seedname)//'.kcontribw',status='unknown')
      open (53,file=trim(seedname)//'.wdoskcontribw',status='unknown')
      open (54,file=trim(seedname)//'.K1w',status='unknown')
      if (resolvecont) then
         open (55,file=trim(seedname)//'.optcondw_resolvedxx',status='unknown')
         open (56,file=trim(seedname)//'.optcondw_resolvedxy',status='unknown')
         open (57,file=trim(seedname)//'.optcondw_resolvedxz',status='unknown')
         open (58,file=trim(seedname)//'.optcondw_resolvedyy',status='unknown')
         open (59,file=trim(seedname)//'.optcondw_resolvedyz',status='unknown')
         open (60,file=trim(seedname)//'.optcondw_resolvedzz',status='unknown')
         allocate(optcond_resolved(nk,1:maxw-1,3,3,nb),totoptcond_resolved(1:maxw-1,3,3,nb))
         allocate(tmpr(nb))
      endif
     
      !read-in old k results
      nkold = 0
       if (.not.band) then
            read(52,*)dum1,nkold,il,jj
       endif
       do jk=1,nkold
         read(52,3000)il,(tmp(jj),jj=1,nsig)
         do jsig=1,nsig
             dccond(jk,sigidx(jsig,1),sigidx(jsig,2)) = tmp(jsig)
         enddo
         read(54,3000)il,(tmp(jj),jj=1,nsig)
         do jsig=1,nsig
             K1(jk,sigidx(jsig,1),sigidx(jsig,2)) = tmp(jsig)
         enddo
         do jw=1,maxw/2-1
           read(52,3001)tmp(1:nsig)
            do jsig=1,nsig
             optcond(jk,jw,sigidx(jsig,1),sigidx(jsig,2)) = tmp(jsig)
           enddo
         enddo
        if (compdos) then
           do jw=-maxw,maxw-1
              read(53,3002)tmpDOS(1:nb)   
              do jb1=1,nb
                 DOS(jk,jw,jb1) = tmpDOS(jb1)
              enddo
           enddo
        endif
      enddo
      close(52)
      close(53)
      close(54)
      

     if (.not.nonint) then
        !   Reading the DMFT Self-Enargy 
        open (553,file='Sfreal.dat',status='old')
        read(553,*)SEnfreq,ninter,intidx(1:ninter)
        write(*,*)SEnfreq,ninter,intidx(1:ninter)
        allocate(SEin(SEnfreq,ninter),wSE(SEnfreq),tmpSE(2*ninter))
        do jw=1,SEnfreq
            read(553,*)wSE(jw),tmpSE!SEin(jw,:)
           do i=1,ninter
             SEin(jw,i) = dcmplx(tmpSE(2*i-1),min(tmpSE(2*i),-deltino))
           enddo
          !write(*,*)wSE(jw),SEin(jw,:)
        enddo
       
        close(553)
        write(*,*)"self energy from input:"
     endif    


      do k=-maxw,maxw-1
         w(k)=real(k)*dw
         do i=1,nb
               SE(k,i) = -Xi*deltino
         enddo
         if (.not.nonint) then
            minimum = 1000d0
            foundw = .false.
            do jw=2,SEnfreq              
               if (.not.(foundw).and.(wSE(jw).gt.w(k))) then
                   foundw = .true.
                   dwSE = wSE(jw)-wSE(jw-1)
                   do i=1,ninter
                      SE(k,intidx(i)) = SEin(jw-1,i)*(wSE(jw)-w(k))/dwSE  &
                                        +SEin(jw,i)*(w(k)-wSE(jw-1))/dwSE
                   enddo
               endif
               !if (dabs(wSE(jw)-w(k)).lt.minimum) then
                   !minimum =  dabs(wSE(jw)-w(k))
                   !do i=1,ninter
                   !   SE(k,intidx(i)) = SEin(jw,i)
                   !enddo
               !endif
            enddo
            if (.not.applyscissors) then
               write(*,"(100F13.9)")w(k),SE(k,:)   
            endif
         endif
         if (applyscissors) then
            do i=1,scissorsnorb
                  SE(k,scissorsorbidx(i)) = SE(k,scissorsorbidx(i)) + dcmplx(scissorsval)
            enddo
            write(*,"(100F13.9)")w(k),SE(k,:)  
         endif
      enddo

      


      diag = .true.
      
      do jk=nkold+1,nk   !loop over k-points
         write(6,*) 'k-point, #=', jk ," of ",nk
         flush(6)

         do jsig=1,nsig  !loop over number of sigma combinations
             idx1 = sigidx(jsig,1)
             idx2 = sigidx(jsig,2)
             do jb1=1,Nb 
                  do jb2=1,Nb 
                       select case (idx1)
                          case (1)
                             matel1(jb1,jb2) = Hkdxfin(jk,jb1,jb2)
                          case (2)
                             matel1(jb1,jb2) = Hkdyfin(jk,jb1,jb2)
                          case (3)
                             matel1(jb1,jb2) = Hkdzfin(jk,jb1,jb2)
                       end select
                       if (idx1.eq.idx2) then 
                           diag = .true.
                       else
                           diag = .false.
                           select case (idx2)
                              case (1)
                                 matel2(jb1,jb2) = Hkdxfin(jk,jb1,jb2)
                              case (2)
                                 matel2(jb1,jb2) = Hkdyfin(jk,jb1,jb2)
                              case (3)
                                 matel2(jb1,jb2) = Hkdzfin(jk,jb1,jb2)
                              end select
                       endif
                   enddo 
             enddo
             do jw=-maxw,maxw-1
                  !set RHS and
                  !obtain LHS=hi=w-Hk-SE 
                  do jb1=1,Nb 
                     do jb2=1,Nb
                        RHS1(jb1,jb2) = matel1(jb1,jb2)
                        RHS1c(jb1,jb2) = matel1(jb1,jb2)
                        LHS(jb1,jb2)=-Hk(jk,jb1,jb2)
                        LHSc(jb1,jb2)=-Hk(jk,jb1,jb2)
                     enddo              
                     LHS(jb1,jb1) = LHS(jb1,jb1) + dcmplx(w(jw)+chempot)-SE(jw,jb1)
                     LHSc(jb1,jb1) = LHSc(jb1,jb1) + dcmplx(w(jw)+chempot)-SE(jw,jb1)
                  enddo 

                  call slse(LHS,RHS1,Nb)
                  call slse_conjg(LHSc,RHS1c,Nb)
                  if ((diag)) then
                     RHS2 = RHS1
                     RHS2c = RHS1c
                  else
                     do jb1=1,Nb 
                        do jb2=1,Nb
                          RHS2(jb1,jb2) = matel2(jb1,jb2)
                          RHS2c(jb1,jb2) = matel2(jb1,jb2)
                        enddo
                     enddo
                     call slse(LHS,RHS2,Nb)
                     call slse_conjg(LHSc,RHS2c,Nb)
                  endif
                  if (compdos.and.(jsig.eq.1)) then
                     do jb1=1,Nb 
                        do jb2=1,Nb
                          eye(jb1,jb2) = 0d0
                        enddo
                        eye(jb1,jb1) = 1d0
                     enddo
                     call slse(LHS,eye,Nb)
                     do jb1=1,Nb 
                        DOS(jk,jw,jb1)=DOS(jk,jw,jb1)-dimag(eye(jb1,jb1))/pi         
                     enddo
                  endif
             
                  do jb1=1,Nb
                     do jb2=1,Nb
                        Gtvk1(jw,jb1,jb2)= (RHS1(jb1,jb2)-RHS1c(jb1,jb2))/2d0/pi/Xi
                        Gtvk2(jw,jb1,jb2)= (RHS2(jb1,jb2)-RHS2c(jb1,jb2))/2d0/pi/Xi

                        enddo
                  enddo

             enddo !frequency loop

             !  computing optical conductivity 
             if(opt) then
                !compute dc conductivity
                do jwint= -maxw,maxw-1
                    if (dabs(beta*w(jwint)).lt.40d0) then
                        df = -beta*dexp(beta*w(jwint))/(dexp(beta*w(jwint))+1d0)**2
                        do jb1=1,Nb
                           do jb2=1,Nb
                              dccond(jk,idx1,idx2) = dccond(jk,idx1,idx2) - dw*df*Gtvk1(jwint,jb1,jb2)*Gtvk2(jwint,jb2,jb1)
                              K1(jk,idx1,idx2) = K1(jk,idx1,idx2) + &
                                             dw*w(jwint)*df*Gtvk1(jwint,jb1,jb2)*Gtvk2(jwint,jb2,jb1)
                              if (resolvecont) then
                                 dccond_resolved(jk,idx1,idx2,jb1) = dccond_resolved(jk,idx1,idx2,jb1) &
                                             - dw*df*Gtvk1(jwint,jb1,jb2)*Gtvk2(jwint,jb2,jb1)
                              K1_resolved(jk,idx1,idx2,jb1) = K1_resolved(jk,idx1,idx2,jb1) + &
                                             dw*w(jwint)*df*Gtvk1(jwint,jb1,jb2)*Gtvk2(jwint,jb2,jb1)
                              endif
                           enddo
                        enddo

                    endif
                enddo

                do jwext= 1,maxw/2-1
                    ! integral over internal frequencies
                    do jwint= -maxw/2,maxw/2
                        
                     !fermi functions
                        if (beta*w(jwint).lt.-40d0) then
                           f_j=1d0  
                        elseif (beta*w(jwint).gt.40d0) then
                           f_j=0d0  
                        else
                           f_j =1.d0/(dexp(beta*w(jwint))+1.d0)
                        endif
                        if (beta*(w(jwint)+w(jwext)).lt.-40d0) then
                           f_kj=1d0  
                        elseif (beta*(w(jwint)+w(jwext)).gt.40d0) then
                           f_kj=0d0  
                        else
                           f_kj =1.d0/(dexp(beta*(w(jwint)+w(jwext)))        +1.d0)
                        endif
                        if (f_j.ne.f_kj) then
				!trace in orbital space
				do jb1=1,Nb
					do jb2=1,Nb
						optcond(jk,jwext,idx1,idx2) = optcond(jk,jwext,idx1,idx2) + dw*(f_j-f_kj)/w(jwext) &
										*Gtvk1(jwint,jb1,jb2)*Gtvk2(jwint+jwext,jb2,jb1)
						if (resolvecont) then
						optcond_resolved(jk,jwext,idx1,idx2,jb1) = optcond_resolved(jk,jwext,idx1,idx2,jb1) &
							+ dw*(f_j-f_kj)/w(jwext)*Gtvk1(jwint,jb1,jb2)*Gtvk2(jwint+jwext,jb2,jb1)
						endif
					enddo
				enddo
                        endif  
                    enddo!internal frequency loop
                enddo!external frequency loop
             endif!opt flag
         enddo!loop over entries of the cond. matrix

!        write(*,*)jk,optcond(jk,1,1,1),Hkdxfin(jk,1,1),Hkdyfin(jk,1,1),Hkdzfin(jk,1,1),Hk(jk,1,1)     
!write(*,*)optcond(jk,1,1,1),optcond(jk,1,1,2),optcond(jk,1,1,3),optcond(jk,1,2,2),optcond(jk,1,2,3),optcond(jk,1,3,3)
      enddo!k-points loop
 
      !write-out results
      if (band) then
         open (52,file=trim(seedname)//'.kcontribw_band',status='unknown')
         write(52,*)"#",nk,maxw/2-1,nsig,convfac
         open (54,file=trim(seedname)//'.K1w_band',status='unknown')
      else
          open (52,file=trim(seedname)//'.kcontribw',status='unknown')
          write(52,*)"#",nk,maxw/2-1,nsig
          open (54,file=trim(seedname)//'.K1w',status='unknown')
      endif
      open (53,file=trim(seedname)//'.wdoskcontribw',status='unknown')
      
      
      do jk=1,nk
         do jsig=1,nsig
             tmp(jsig) = dccond(jk,sigidx(jsig,1),sigidx(jsig,2))
         enddo
         write(52,3000)jk,(tmp(jj),jj=1,nsig)
         do jsig=1,nsig
             tmp(jsig) = K1(jk,sigidx(jsig,1),sigidx(jsig,2))
         enddo
         write(54,3000)jk,(tmp(jj),jj=1,nsig)
         do jw=1,maxw/2-1
           do jsig=1,nsig
             tmp(jsig) = optcond(jk,jw,sigidx(jsig,1),sigidx(jsig,2))
           enddo
           write(52,3001)tmp(1:nsig)
         enddo
         if (compdos) then
           do jw=-maxw,maxw-1
              do jb1=1,nb
                 tmpDOS(jb1) = DOS(jk,jw,jb1)
              enddo
              write(53,3002)tmpDOS(1:nb)
           enddo
        endif
      enddo
      close(52)
      close(53)
      close(54)
      
3000 FORMAT(I8,6E20.12)
3001 FORMAT(6E20.12)
3002 FORMAT(100E20.12)

     if (band) then
         stop
     endif

     do jk=1,nk
         do jd1=1,3
             do jd2=jd1+1,3
                  dccond(jk,jd2,jd1) = dccond(jk,jd1,jd2)
                  K1(jk,jd2,jd1) = K1(jk,jd1,jd2)
                  do jw=1,maxw/2-1
                     optcond(jk,jw,jd2,jd1) = optcond(jk,jw,jd1,jd2)
                  enddo
                  if (resolvecont) then
                      do jb1=1,nb
			dccond_resolved(jk,jd2,jd1,jb1) = dccond_resolved(jk,jd1,jd2,jb1)
			K1_resolved(jk,jd2,jd1,jb1) = K1_resolved(jk,jd1,jd2,jb1)
			do jw=1,maxw/2-1
			  optcond_resolved(jk,jw,jd2,jd1,jb1) = optcond_resolved(jk,jw,jd1,jd2,jb1)
  		        enddo
                      enddo
                  endif
             enddo
         enddo
     enddo
     totoptcond = 0d0
     totdccond = 0d0
     totK1 = 0d0
     if (resolvecont) then
         totoptcond_resolved = 0d0
         totdccond_resolved = 0d0
         totK1_resolved = 0d0
     endif
     
     open(unit=44,file=trim(seedname)//".symop")
     read(44,*)nop
     allocate(symop(nop,3,3))
     do js=1,nop
        do jd=1,3
           read(44,"(3F8.5)")symop(js,jd,:)
!            write(*,*)js,lattice%sym(js,jd,:)
        enddo
        read(44,*)
     enddo
     close (44)

      ! symmetrized tetrahedral integration
      write(*,*) ">>> symmetrized tetrahedral integration"
      sumtetraweights = 0d0
      do jt=1,nt
          sumtetraweights = sumtetraweights + wtetra(jt)
      enddo     
      do jt=1,nt
         reftet(1,1:4) = (/tetra(jt,1),tetra(jt,5),tetra(jt,6),tetra(jt,7)/)
	 reftet(2,1:4) = (/tetra(jt,2),tetra(jt,5),tetra(jt,8),tetra(jt,9)/)
	 reftet(3,1:4) = (/tetra(jt,3),tetra(jt,6),tetra(jt,8),tetra(jt,10)/)
	 reftet(4,1:4) = (/tetra(jt,4),tetra(jt,7),tetra(jt,9),tetra(jt,10)/)
! 	 if (tclass(jt).eq.1) then !central octahedron
		  !corresp. to class 1 in Endres paper
	 reftet(5,1:4) = (/tetra(jt,5),tetra(jt,6),tetra(jt,7),tetra(jt,10)/)
	 reftet(6,1:4) = (/tetra(jt,5),tetra(jt,6),tetra(jt,8),tetra(jt,10)/)
	 reftet(7,1:4) = (/tetra(jt,5),tetra(jt,7),tetra(jt,9),tetra(jt,10)/)
	 reftet(8,1:4) = (/tetra(jt,5),tetra(jt,8),tetra(jt,9),tetra(jt,10)/)
! 	 else
! 	      !corresp. to class 2 in Endres paper
! 	      reftet(5,1:4) = (/tetra(jt,5),tetra(jt,6),tetra(jt,7),tetra(jt,9)/)
! 	      reftet(6,1:4) = (/tetra(jt,5),tetra(jt,6),tetra(jt,8),tetra(jt,9)/)
! 	      reftet(7,1:4) = (/tetra(jt,6),tetra(jt,7),tetra(jt,9),tetra(jt,10)/)
! 	      reftet(8,1:4) = (/tetra(jt,6),tetra(jt,8),tetra(jt,9),tetra(jt,10)/)
! 	 endif 
         do js=1,nop
             do jd1=1,3
                 do jd2=1,3
                    do jd3=1,3
                       do jd4=1,3
!                            do jd=1,4
!                               totdccond(jd1,jd2) = totdccond(jd1,jd2) &
!                      - symop(js,jd3,jd1)*dccond(tetra(jt,jd),jd3,jd4)*symop(js,jd4,jd2)*wtetra(jt)/20d0/dble(nop)
!                                totK1(jd1,jd2) = totK1(jd1,jd2) &
!                      - symop(js,jd3,jd1)*K1(tetra(jt,jd),jd3,jd4)*symop(js,jd4,jd2)*wtetra(jt)/20d0/dble(nop)
!                               do jw=1,maxw/2-1
!                                     totoptcond(jw,jd1,jd2) = totoptcond(jw,jd1,jd2) &
!                      - symop(js,jd3,jd1)*optcond(tetra(jt,jd),jw,jd3,jd4)*symop(js,jd4,jd2)*wtetra(jt)/20d0/dble(nop)
!                               enddo
!                            enddo
!                            do jd=5,10
!                               totdccond(jd1,jd2) = totdccond(jd1,jd2) &
!                      + symop(js,jd3,jd1)*dccond(tetra(jt,jd),jd3,jd4)*symop(js,jd4,jd2)*wtetra(jt)/5d0/dble(nop)
!                               totK1(jd1,jd2) = totK1(jd1,jd2) &
!                      + symop(js,jd3,jd1)*K1(tetra(jt,jd),jd3,jd4)*symop(js,jd4,jd2)*wtetra(jt)/5d0/dble(nop)
!                               do jw=1,maxw/2-1
!                                     totoptcond(jw,jd1,jd2) = totoptcond(jw,jd1,jd2) &
!                      + symop(js,jd3,jd1)*optcond(tetra(jt,jd),jw,jd3,jd4)*symop(js,jd4,jd2)*wtetra(jt)/5d0/dble(nop)
!                               enddo
!                            enddo
                          do jd=1,8
                           do jda=1,4
                              totdccond(jd1,jd2) = totdccond(jd1,jd2) &
                     + symop(js,jd3,jd1)*dccond(reftet(jd,jda),jd3,jd4)*symop(js,jd4,jd2)*wtetra(jt)/32d0/dble(nop)
                              totK1(jd1,jd2) = totK1(jd1,jd2) &
                     + symop(js,jd3,jd1)*K1(reftet(jd,jda),jd3,jd4)*symop(js,jd4,jd2)*wtetra(jt)/32d0/dble(nop)
                              do jw=1,maxw/2-1
                                    totoptcond(jw,jd1,jd2) = totoptcond(jw,jd1,jd2) &
                     + symop(js,jd3,jd1)*optcond(reftet(jd,jda),jw,jd3,jd4)*symop(js,jd4,jd2)*wtetra(jt)/32d0/dble(nop)
                              enddo
                           enddo
                         enddo
                           if (resolvecont) then 
                                do jb1=1,nb
				    do jd=1,4
					  totdccond_resolved(jd1,jd2,jb1) = totdccond_resolved(jd1,jd2,jb1) &
				- symop(js,jd3,jd1)*dccond_resolved(tetra(jt,jd),jd3,jd4,jb1)*symop(js,jd4,jd2)*wtetra(jt)/20d0/dble(nop)
					  totK1_resolved(jd1,jd2,jb1) = totK1_resolved(jd1,jd2,jb1) &
				- symop(js,jd3,jd1)*K1_resolved(tetra(jt,jd),jd3,jd4,jb1)*symop(js,jd4,jd2)*wtetra(jt)/20d0/dble(nop)
					  do jw=1,maxw/2-1
						totoptcond_resolved(jw,jd1,jd2,jb1) = totoptcond_resolved(jw,jd1,jd2,jb1) &
				- symop(js,jd3,jd1)*optcond_resolved(tetra(jt,jd),jw,jd3,jd4,jb1)*symop(js,jd4,jd2)*wtetra(jt)/20d0/dble(nop)
					  enddo
				      enddo
				      do jd=5,10
					  totdccond_resolved(jd1,jd2,jb1) = totdccond_resolved(jd1,jd2,jb1) &
				+ symop(js,jd3,jd1)*dccond_resolved(tetra(jt,jd),jd3,jd4,jb1)*symop(js,jd4,jd2)*wtetra(jt)/5d0/dble(nop)
					  totK1_resolved(jd1,jd2,jb1) = totK1_resolved(jd1,jd2,jb1) &
				+ symop(js,jd3,jd1)*K1_resolved(tetra(jt,jd),jd3,jd4,jb1)*symop(js,jd4,jd2)*wtetra(jt)/5d0/dble(nop)
					  do jw=1,maxw/2-1
						totoptcond_resolved(jw,jd1,jd2,jb1) = totoptcond_resolved(jw,jd1,jd2,jb1) &
				+ symop(js,jd3,jd1)*optcond_resolved(tetra(jt,jd),jw,jd3,jd4,jb1)*symop(js,jd4,jd2)*wtetra(jt)/5d0/dble(nop)
					  enddo
				      enddo
                                enddo      
                           endif
                       enddo
                    enddo
                 enddo
             enddo
            !write(*,*)js,totdccond(1,1),totoptcond(1,1,1)
         enddo
         if (compdos) then
!              do jd1=1,nb
!                 do jw=-maxw,maxw-1
!                     do jd=1,4
!                        DOStot(jw,jd1) = DOStot(jw,jd1) - DOS(tetra(jt,jd),jw,jd1)*wtetra(jt)/20d0
!                     enddo
!                     do jd=5,10
!                        DOStot(jw,jd1) = DOStot(jw,jd1) + DOS(tetra(jt,jd),jw,jd1)*wtetra(jt)/5d0
!                     enddo
!                 enddo
!              enddo
              do jd1=1,nb
                do jw=-maxw,maxw-1
                    do jd=1,8
                           do jda=1,4
                               DOStot(jw,jd1) = DOStot(jw,jd1) + DOS(reftet(jd,jda),jw,jd1)*wtetra(jt)/32d0
                           enddo                
                    enddo
                enddo
             enddo
         endif
      enddo

      !compute thermopower
      thermopower = 0d0
      !explanation of factor: conv2muV*kboltz/echarge
      do jd1=1,3
         do jd2=1,3
             thermopower(jd1,jd2) = 86.173427909006634d0*totK1(jd1,jd2)/totdccond(jd1,jd2)*beta
             if (resolvecont) then
                 do jb1=1,nb
                     thermopower_resolved(jd1,jd2,jb1) = &
                         86.173427909006634d0*totK1_resolved(jd1,jd2,jb1)/totdccond_resolved(jd1,jd2,jb1)*beta
                 enddo
             endif
         enddo
      enddo
      

      !compute DOS
      if (compdos) then   
         do jd1=1,nb
            do jw=-maxw,maxw-1
                totDOS(jw) = totDOS(jw) + DOStot(jw,jd1)
            enddo
         enddo
      endif

      !unit conversion to Siemens/cm
      do jd1=1,3
          do jd2=1,3
             totdccond(jd1,jd2) = convfac*totdccond(jd1,jd2)
             do jw=1,maxw/2-1
                totoptcond(jw,jd1,jd2) = convfac*totoptcond(jw,jd1,jd2)
             enddo
             if (resolvecont) then
                 do jb1=1,nb
                    totdccond_resolved(jd1,jd2,jb1) = convfac*totdccond_resolved(jd1,jd2,jb1)
		    do jw=1,maxw/2-1
			totoptcond_resolved(jw,jd1,jd2,jb1) = convfac*totoptcond_resolved(jw,jd1,jd2,jb1)
		    enddo
                 enddo 
             endif
          enddo
      enddo
      write(6,*)'dc cond. xx,yy,zz [Ohm-1 cm-1]',totdccond(1,1),totdccond(2,2),totdccond(3,3)
      write(6,*)'thermopower xx,yy,zz [muV/K]',thermopower(1,1),thermopower(2,2),thermopower(3,3)
      
      write(6,*) 'Step for intergral=', dw
      write(6,*) 'sum of tetrahedral weights=', sumtetraweights
      sumrule = 0d0
      sumruledrude = 0d0
      do jsig=1,nsig
          tmp(jsig) = totdccond(sigidx(jsig,1),sigidx(jsig,2))
      enddo
      !write header for optcondw
      write(50,*)"# Optical conductivity in Ohm-1 cm-1"
      write(50,*)"# w  sigma_xx sigma_xy sigma_xz sigma_yy sigma_yz sigma_zz"
      write(50,'(1E20.12,100E20.12)') 0d0,(tmp(jj),jj=1,nsig) 
      do jw=1,maxw/2-1
         counter = 0
         do jsig=1,nsig
             if (w(jw).le.drudesep) then
                sumruledrude(jsig) = sumruledrude(jsig) + totoptcond(jw,sigidx(jsig,1),sigidx(jsig,2))*dw
             else
                sumrule(jsig) = sumrule(jsig) + totoptcond(jw,sigidx(jsig,1),sigidx(jsig,2))*dw
             endif
         enddo
         do jsig=1,nsig
             tmp(jsig) = totoptcond(jw,sigidx(jsig,1),sigidx(jsig,2))
         enddo
         write(50,'(1E20.12,100E20.12)') w(jw),(tmp(jj), jj=1,nsig)
      enddo
      write(6,*)"conversion factor:",convfac
      write(6,*)"sumrules1[Drude,interband]:",sumruledrude(1),sumrule(1)
      do jsig=1,nsig       
         sumruledrude(jsig) = sumruledrude(jsig) + totdccond(sigidx(jsig,1),sigidx(jsig,2))/sumtetraweights*dw
      enddo
      write(6,*)"sumrules2[Drude,interband]:",sumruledrude(1),sumrule(1)
      close(50)
      if (compdos) then
        write(51,*)"# spectral function in eV^-1"
        write(51,*)"# w  total orbital1 orbital2 ...."
        do jw=-maxw,maxw-1
            write(51,'(200F13.8)')w(jw),totDOS(jw),DOStot(jw,:)
        enddo
      endif
      close(51)
      close(52)
     
      sumruledrude1 = sumruledrude(1)
      sumrule1 = sumrule(1)
      
      if (resolvecont) then
         
         do jsig=1,nsig
             do jb1=1,nb
                tmpr(jb1) = thermopower_resolved(sigidx(jsig,1),sigidx(jsig,2),jb1)
             enddo
             write(54+jsig,'(A2,100E20.12)') "# ",(tmpr(jj), jj=1,nb)
             do jb1=1,nb
                tmpr(jb1) = totdccond_resolved(sigidx(jsig,1),sigidx(jsig,2),jb1)
             enddo
             write(54+jsig,'(1E20.12,100E20.12)')0d0,(tmpr(jj), jj=1,nb)
             do jw=1,maxw/2-1
                  do jb1=1,nb
		      tmpr(jb1) = totoptcond_resolved(jw,sigidx(jsig,1),sigidx(jsig,2),jb1)
		  enddo
		  write(54+jsig,'(1E20.12,100E20.12)') w(jw),(tmpr(jj), jj=1,nb)
             enddo
             close(54+jsig)
         enddo
      endif
      return
      end
      
      
! c     ****************************************************************** 
SUBROUTINE slse(A,B,N)
 !solves a system of linear equations with lapack routines
 !A is the LHS matrix and B is the matrix of RHS vectors,
 !B is overwritten and contains the solution vectors on output
 implicit none
 
 integer :: N ,INFO
 integer :: ipiv(N),j1,j2
 complex*16 :: A(N,N),B(N,N),tmp(N,N)
 
 do j1=1,N
   do j2=1,N
      tmp(j1,j2) = A(j1,j2)
   enddo
 enddo
 call zgetrf(N,N,tmp,N,ipiv,INFO)
 if (INFO.ne.0) then
   write(*,*)"ERROR in lapack zgetrf, INFO=",INFO
   stop
 endif
 call zgetrs('N',N,N,tmp,N,ipiv,B,N,INFO)
 if (INFO.ne.0) then
   write(*,*)"ERROR in lapack zgetrs, INFO=",INFO
   stop
 endif


 END SUBROUTINE

! c     ****************************************************************** 
SUBROUTINE slse_conjg(A,B,N)
 !solves a system of linear equations with lapack routines
 !A is the LHS matrix and B is the matrix of RHS vectors,
 !B is overwritten and contains the solution vectors on output
 !conjugated version
 implicit none
 
 integer :: N ,INFO
 integer :: ipiv(N),j1,j2
 complex*16 :: A(N,N),B(N,N),tmp(N,N)
 
 do j1=1,N
   do j2=1,N
      tmp(j1,j2) = A(j1,j2)
   enddo
 enddo
 call zgetrf(N,N,tmp,N,ipiv,INFO)
 if (INFO.ne.0) then
   write(*,*)"ERROR in lapack zgetrf, INFO=",INFO
   stop
 endif
 call zgetrs('C',N,N,tmp,N,ipiv,B,N,INFO)
 if (INFO.ne.0) then
   write(*,*)"ERROR in lapack zgetrs, INFO=",INFO
   stop
 endif


 END SUBROUTINE

   integer function line_count(fid)  
   !>Returns the number of lines in a file.
   !>\param fname Name of the file
   !>\return Number of lines in the file.
	
   implicit none

   !input parameters
!    character(len=*) fname   !filename of the file to count
   integer fid 

   !local parameters
   character*20 dummy
   integer ioStat
   logical ioEndFlag

   ioEndFlag = .false.
   line_count = 0
   do while (.not. ioEndFlag )
      read(fid,"(A20)", iostat=ioStat ) dummy
      if( iostat .eq. 0) line_count = line_count + 1
      if( iostat < 0 ) ioEndFlag = .TRUE. 
   end do
   rewind(fid)

    END FUNCTION line_count   


  SUBROUTINE get_umatrix(H,nw,U)

  implicit none
 
  integer nw,info,j,i
  complex*16 H(nw,nw),work(2*nw),u(nw,nw)
  real*8 eigval(nw),rwork(7*nw)

   do i=1,nw
       do j=1,nw
           U(i,j) = H(i,j)
       enddo
   enddo
  call ZHEEV( 'V', 'U', nw, U, nw, eigval, work, 2*nw, rwork,info )

  if (info.ne.0) then
    write(*,*)"Error: in zheev lapack routine, info:",info
    do i=1,nw
       do j=1,nw
          write(*,*)i,j,H(i,j)
       enddo
    enddo
    stop
  endif


  END SUBROUTINE

 
 SUBROUTINE nonuniform_stencil(points,dist,outp)
 !computes the derivative along an axis with 5 non-uniform points

 implicit none

 real*8 ::  dist(4),c1,c2,c3
 complex*16 :: points(5),outp

   c1 = dist(3)**2*(dist(3)+dist(2))&
               /(dist(3)+dist(4))**2 &           
               /(dist(1)+dist(2)+dist(3)+dist(4))
   c2 = 1d0-dist(3)**2/dist(2)**2-c1+&
                c1*(dist(3)+dist(4))**2/(dist(2)+dist(1))**2
   c3 = dist(3)+dist(3)**2/dist(2)-&
             c1*(dist(3)+dist(4))*(1d0+(dist(3)+dist(4))/(dist(2)+dist(1)))
    
    outp = points(4)- cmplx(dist(3)**2/dist(2)**2)*points(2)&
                    - cmplx(c1)*points(5) &
                    + cmplx(c1*((dist(3)+dist(4))/(dist(2)+dist(1)))**2)*points(1)&
                    - cmplx(c2)*points(3)

    outp = outp/c3;

 END SUBROUTINE

 SUBROUTINE debugfunct3(seedname,nk,k,nt,tetra,wtetra,lattice,ndim)
!fermi derivative at (kx,ky,kz)=0  f(x)=-exp(kx)/(exp(kx)+1)**2 f_tot=f(x)f(y)f(z)
!result :-0.014691482994426 (=(2/(e^(1/2)+1)-1/2)^3 for beta=1d0 (cubic)
       ! -0.960378027086140 for beta 10d0(cubic)
       ! -0.999999987633078 for beta 40d0(cubic)
       ! -1 for beta 100d0(cubic)
       use util 

      implicit none

      integer nk,nt,i,jd,jk,j1,j2,ndim(3)
      integer k(nk,3),tetra(nt,10),ktmp(nk,3),reftet(8,4)
      real*8 wtetra(nt),kcontrib(nk),sumtetraweights,res,maxk,loccont,res2,res3
      real*8 loccont2,beta
      real*8 rotmat(3,3),rvec(3)
      type(structure) lattice
 
      character*50 seedname

!        reftet(1,1:4) = (/tetra(1),tetra(5),tetra(6),tetra(7)/)
!        reftet(2,1:4) = (/tetra(2),tetra(5),tetra(8),tetra(9)/)
!        reftet(3,1:4) = (/tetra(3),tetra(6),tetra(8),tetra(10)/)
!        reftet(4,1:4) = (/tetra(4),tetra(7),tetra(9),tetra(10)/)
!        reftet(5,1:4) = (/tetra(5),tetra(6),tetra(7),tetra(10)/)
!        reftet(6,1:4) = (/tetra(5),tetra(6),tetra(8),tetra(10)/)
!        reftet(7,1:4) = (/tetra(5),tetra(7),tetra(9),tetra(10)/)
!        reftet(8,1:4) = (/tetra(5),tetra(8),tetra(9),tetra(10)/)

      call inverse3x3(lattice%transform_matrix,rotmat)
      ktmp = 0
        write(*,*)"ndim",ndim
       write(*,*)"1",rotmat(1,:)
       write(*,*)"2",rotmat(2,:)
       write(*,*)"3",rotmat(3,:)
!       stop
      do jk=1,nk
         do j1=1,3
            rvec(j1) = real(k(jk,j1))/real(ndim(j1))
         enddo
         if ((lattice%specifier.ne."R").and.(lattice%specifier.ne."H")) then
            rvec = matmul(lattice%transform_matrix,rvec)
         endif
         do j1=1,3
            ktmp(jk,j1) = int(rvec(j1)*ndim(j1))
         enddo
         write(*,*)"kj",jk,k(jk,:),ktmp(jk,:)
      enddo

      k = ktmp

      beta = 40d0
      sumtetraweights = 0d0
      kcontrib = 0d0
      res = 0d0
      res2 = 0d0
      res3 = 0d0
      maxk = 0d0
!       do jk=1,nk
!          if (dble(k(jk,1)).gt.maxk) then
             maxk = ndim(1)!dble(k(jk,1))
!          endif
!       enddo
      do jk=1,nk
!           kcontrib(jk)
!          do jd=1,3
!            if (real(k(jk,jd))/maxk.lt.2d0/3d0) then
!               kcontrib(jk) = kcontrib(jk) + (real(k(jk,jd))/maxk)**2*9d0/4d0
!            elseif (real(k(jk,jd))/maxk.gt.2d0/3d0) then
!               kcontrib(jk) = kcontrib(jk) + ((maxk - real(k(jk,jd)))/maxk)**2*9d0
!            else
!               kcontrib(jk) = kcontrib(jk) + 1d0
!             endif
            kcontrib(jk) = -beta*dexp(beta*dble(k(jk,1))/maxk)/(dexp(beta*dble(k(jk,1))/maxk)+1d0)**2*&
                            beta*dexp(beta*dble(k(jk,2))/maxk)/(dexp(beta*dble(k(jk,2))/maxk)+1d0)**2*&
                            beta*dexp(beta*dble(k(jk,3))/maxk)/(dexp(beta*dble(k(jk,3))/maxk)+1d0)**2
!          enddo
         write(*,*)k(jk,:),kcontrib(jk)
      enddo

      do i=1,nt
           sumtetraweights = sumtetraweights + wtetra(i)
      enddo
      do i=1,nt
         reftet(1,1:4) = (/tetra(i,1),tetra(i,5),tetra(i,6),tetra(i,7)/)
         reftet(2,1:4) = (/tetra(i,2),tetra(i,5),tetra(i,8),tetra(i,9)/)
         reftet(3,1:4) = (/tetra(i,3),tetra(i,6),tetra(i,8),tetra(i,10)/)
         reftet(4,1:4) = (/tetra(i,4),tetra(i,7),tetra(i,9),tetra(i,10)/)
         reftet(5,1:4) = (/tetra(i,5),tetra(i,6),tetra(i,7),tetra(i,10)/)
         reftet(6,1:4) = (/tetra(i,5),tetra(i,6),tetra(i,8),tetra(i,10)/)
         reftet(7,1:4) = (/tetra(i,5),tetra(i,7),tetra(i,9),tetra(i,10)/)
         reftet(8,1:4) = (/tetra(i,5),tetra(i,8),tetra(i,9),tetra(i,10)/)
        loccont = 0d0
        loccont2 = 0d0
         do jd=1,4
            res = res - kcontrib(tetra(i,jd))*wtetra(i)/20d0
           loccont = loccont- kcontrib(tetra(i,jd))*wtetra(i)/20d0
           loccont2 = loccont2+ kcontrib(tetra(i,jd))*wtetra(i)/4d0
           res2 = res2 + kcontrib(tetra(i,jd))*wtetra(i)/4d0
!            res3 = res3 + kcontrib(tetra(i,jd))*wtetra(i)/32d0
         enddo
         do jd=5,10
            res = res + kcontrib(tetra(i,jd))*wtetra(i)/5d0
             loccont = loccont+ kcontrib(tetra(i,jd))*wtetra(i)/5d0
         enddo
         do j1=1,8
            do jd=1,4
               res3 = res3 + kcontrib(reftet(j1,jd))*wtetra(i)/32d0
            enddo
         enddo
         write(*,*)wtetra(i),tetra(i,:),loccont,loccont2
      enddo

      write(6,*)"sumoweights: ",sumtetraweights
      write(6,*)"res: ",res,res2,res3
!write-out results
      open (52,file=trim(seedname)//'.kcontribw',status='unknown')
      open (53,file=trim(seedname)//'.wdoskcontribw',status='unknown')
      open (54,file=trim(seedname)//'.K1w',status='unknown')
      write(52,*)"#",nk,0,1
      do jk=1,nk
!          do jsig=1,nsig
!              tmp(jsig) = kcontrib(jk)
!          enddo
         write(52,3000)jk,kcontrib(jk)
!          do jsig=1,nsig
!              tmp(jsig) = K1(jk,sigidx(jsig,1),sigidx(jsig,2))
!          enddo
!          write(54,3000)jk,(tmp(jj),jj=1,nsig)
!          do jw=1,maxw/2-1
!            do jsig=1,nsig
!              tmp(jsig) = optcond(jk,jw,sigidx(jsig,1),sigidx(jsig,2))
!            enddo
!            write(52,3001)tmp(1:nsig)
!          enddo
!          if (compdos) then
!            do jw=-maxw,maxw-1
!               do jb1=1,nb
!                  tmpDOS(jb1) = DOS(jk,jw,jb1)
!               enddo
!               write(53,3002)tmpDOS(1:nb)
!            enddo
!         endif
      enddo
      close(52)
      close(53)
      close(54)
      return
3000 FORMAT(I8,6E20.12)
      end

 SUBROUTINE compute_tetrahedralderivative(nk,tetra,k,dk,Hloctet,node,kdim,derivative,debmode,rbasis)
 !computes the derivative within a tetrahedron along cartesian coordinate dimidx

 implicit none
  
 integer nk,k(nk,3),dimidx,tetra(10),dNdt(10,4,10)
 integer j1,j2,j3,ipiv(4),node,INFO,intmap(10),debmode
 integer jn1,jn2,jd1,kdim(3),jd2
 complex*16 :: Hloctet(10),derivative(3),RHS(4,3),func(4),A(4,4)
 real*8 :: dk(3),lock(10,3),rbasis(3,3),tmp(3)

 do jn1=1,10
   lock(jn1,:) = dble(k(tetra(jn1),:))
 enddo
 do jn1=1,10
    do jn2=1,10
       do jd1=1,3
          if ((k(tetra(jn1),jd1).eq.0).and.(k(tetra(jn2),jd1).gt.kdim(jd1)/2)) then
             lock(jn1,jd1) = dble(kdim(jd1))
             write(*,*)"warning:there might be something about the k-basis"
          endif
       enddo
    enddo
 enddo 
 if (debmode.eq.1) then
      do j1=1,10
      write(*,*)"bef",tetra(j1),lock(j1,:)
       enddo
   endif

  do jn1=1,10
     tmp = lock(jn1,:)
     lock(jn1,:) = 0d0
     do jd1=1,3
         do jd2=1,3
            lock(jn1,jd1) = lock(jn1,jd1) + rbasis(jd1,jd2)*tmp(jd2)/dble(kdim(jd2))
         enddo
     enddo
  enddo

if (debmode.eq.1) then
 write(*,*)"this is your input",nk,dk,node,kdim
 do j1=1,10
   write(*,*)"aft",tetra(j1),lock(j1,:)
 enddo
endif
if (debmode.eq.1) then
 do j1=1,10
   write(*,*)"Hloc",tetra(j1),Hloctet(j1)
 enddo
endif
!assemble dNdt
 dNdt = 0
!at node1
dNdt(1,1,1)=3;
dNdt(2,2,1)=-1;
dNdt(5,2,1)=4;
dNdt(3,3,1)=-1;
dNdt(7,3,1)=4;
dNdt(4,4,1)=-1;
dNdt(8,4,1)=4;
!at node2
dNdt(1,1,2)=-1;
dNdt(5,1,2)=4;
dNdt(2,2,2)=3;
dNdt(3,3,2)=-1;
dNdt(6,3,2)=4;
dNdt(4,4,2)=-1;
dNdt(9,4,2)=4;
!at node3
dNdt(1,1,3)=-1;
dNdt(7,1,3)=4;
dNdt(2,2,3)=-1;
dNdt(6,2,3)=4;
dNdt(3,3,3)=3;
dNdt(4,4,3)=-1;
dNdt(10,4,3)=4;
!at node4
dNdt(1,1,4)=-1;
dNdt(8,1,4)=4;
dNdt(2,2,4)=-1;
dNdt(9,2,4)=4;
dNdt(3,3,4)=-1;
dNdt(10,3,4)=4;
dNdt(4,4,4)=3;
!at node5
dNdt(5,1,5)=2;
dNdt(1,1,5)=1;
dNdt(2,2,5)=1;
dNdt(5,2,5)=2;
dNdt(3,3,5)=-1;
dNdt(6,3,5)=2;
dNdt(7,3,5)=2;
dNdt(4,4,5)=-1;
dNdt(8,4,5)=2;
dNdt(9,4,5)=2;
!at node6
dNdt(1,1,6)=-1;
dNdt(5,1,6)=2;
dNdt(7,1,6)=2;
dNdt(2,2,6)=1;
dNdt(6,2,6)=2;
dNdt(3,3,6)=1;
dNdt(6,3,6)=2;
dNdt(4,4,6)=-1;
dNdt(9,4,6)=2;
dNdt(10,4,6)=2;
!at node7
dNdt(1,1,7)=1;
dNdt(7,1,7)=2;
dNdt(2,2,7)=-1;
dNdt(5,2,7)=2;
dNdt(6,2,7)=2;
dNdt(3,3,7)=1;
dNdt(7,3,7)=2;
dNdt(4,4,7)=-1;
dNdt(8,4,7)=2;
dNdt(10,4,7)=2;
!at node8
dNdt(1,1,8)=1;
dNdt(8,1,8)=2;
dNdt(2,2,8)=-1;
dNdt(5,2,8)=2;
dNdt(9,2,8)=2;
dNdt(3,3,8)=-1;
dNdt(7,3,8)=2;
dNdt(10,3,8)=2;
dNdt(4,4,8)=1;
dNdt(8,4,8)=2;
!at node9
dNdt(1,1,9)=-1;
dNdt(5,1,9)=2;
dNdt(8,1,9)=2;
dNdt(2,2,9)=1;
dNdt(9,2,9)=2;
dNdt(3,3,9)=-1;
dNdt(6,3,9)=2;
dNdt(10,3,9)=2;
dNdt(4,4,9)=1;
dNdt(9,4,9)=2;
!at node10
dNdt(1,1,10)=-1;
dNdt(7,1,10)=2;
dNdt(8,1,10)=2;
dNdt(2,2,10)=-1;
dNdt(6,2,10)=2;
dNdt(9,2,10)=2;
dNdt(3,3,10)=1;
dNdt(10,3,10)=2;
dNdt(4,4,10)=1;
dNdt(10,4,10)=2;

intmap = (/1,2,3,4,5,7,8,6,9,10/)
func = cmplx(0d0)
do j1=1,4
  do j2=1,10
    func(j1) = func(j1) + cmplx(dNdt(intmap(j2),j1,intmap(node)))*Hloctet(j2)
  enddo
 if (debmode.eq.1) then
     write(*,*)"func",j1,func(j1)
  endif
enddo
 
A = 0
A(1,:) = cmplx(1d0)
do j1=1,3
  do j2=1,4
     do j3=1,10
       if (debmode.eq.1) then
       endif
       A(1+j1,j2) = A(1+j1,j2) + cmplx(dNdt(intmap(j3),j2,intmap(node)))*cmplx(lock(j3,j1))!*cmplx(dk(j1))
     enddo
  enddo
  if (debmode.eq.1) then

  endif
enddo
RHS = cmplx(0d0)
do j1=1,3
  RHS(j1+1,j1) = cmplx(1d0)
enddo

call zgetrf(4,4,A,4,ipiv,INFO)
if (INFO.ne.0) then
  write(*,*)"ERROR in lapack zgetrf, INFO=",INFO
  write(*,*)"A is,",A(1,:)
  write(*,*)"A is,",A(2,:)
  write(*,*)"A is,",A(3,:)
  write(*,*)"A is,",A(4,:)
  write(*,*)"k is,",lock(1,:)
  write(*,*)"k is,",lock(2,:)
  write(*,*)"k is,",lock(3,:)
  write(*,*)"k is,",lock(4,:)
  write(*,*)"k is,",lock(5,:)
  write(*,*)"k is,",lock(6,:)
  write(*,*)"k is,",lock(7,:)
  write(*,*)"k is,",lock(8,:)
  write(*,*)"k is,",lock(9,:)
  write(*,*)"k is,",lock(10,:)
  write(*,*)"ko is,",k(tetra(1),:)
  write(*,*)"ko is,",k(tetra(2),:)
  write(*,*)"ko is,",k(tetra(3),:)
  write(*,*)"ko is,",k(tetra(4),:)
  write(*,*)"ko is,",k(tetra(5),:)
  write(*,*)"ko is,",k(tetra(6),:)
  write(*,*)"ko is,",k(tetra(7),:)
  write(*,*)"ko is,",k(tetra(8),:)
  write(*,*)"ko is,",k(tetra(9),:)
  write(*,*)"ko is,",k(tetra(10),:)
  write(*,*)"tetra is",tetra
  write(*,*)"rbasis is",rbasis(1,:)
  write(*,*)"rbasis is",rbasis(2,:)
  write(*,*)"rbasis is",rbasis(3,:)
  stop
endif
call zgetrs('N',4,3,A,4,ipiv,RHS,4,INFO)
if (INFO.ne.0) then
  write(*,*)"ERROR in lapack zgetrs, INFO=",INFO
  stop
endif
if (debmode.eq.1) then
  do j1=1,4
    write(*,*)"RHS",RHS(j1,:)
  enddo
endif
derivative = cmplx(0d0)

do j1=1,4
   do j2=1,3
      derivative(j2) = derivative(j2) + RHS(j1,j2)*func(j1)
   enddo
enddo
if (debmode.eq.1) then
   write(*,*)"der",derivative
endif

 END SUBROUTINE

