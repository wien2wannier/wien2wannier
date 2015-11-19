  PROGRAM refine_tmesh
  !adaptive refinement of tetrahedral k-mesh 

  use util 

 implicit none

 integer :: nk,jk,nkfull,jd,find_kpoint2,lock(3),idiv,ndiv(3),nt,locv
 integer, allocatable :: k(:,:),kfull(:,:),kfulltmp(:,:)!,nnk(:,:),kdist(:,:)
 integer :: maxidx,kmax(3),jkfull,jv,tchange
 integer,allocatable :: map(:,:),tclass(:)
 integer :: rsteps,counter,maximum,nnewk,ir,nnewweights,initsteps,roundi
 real*8 :: theta,sumw
 real*8, allocatable :: kcontribw(:,:,:),nabk(:),tmp(:),kfull_int(:,:)
 integer :: kidx(6),jd1,jd2,npidx(12,2),nev,eidx1,eidx2,eidx11,eidx22
 integer, allocatable :: tetra(:,:)
 real*8, allocatable :: wtetra(:)
 integer, allocatable :: newk(:,:),newtetra(:,:),newtclass(:)
 integer, allocatable :: VOE(:,:,:),VOEidx(:,:)
 real*8 :: locmean,maxvariance
 real*8 :: rdum
 integer :: nnewt,tcount,newpoints(6,3),nop,nvoe
 logical,allocatable :: trefine(:)
 real*8, allocatable :: newwtetra(:),tvariance(:)
 integer :: js,loct(4),nsig,nw,jw
 character*80 dummy
 character*80 seedname
 character*80 argdummy
 integer IK,kdiv,nkp,kdim(3),i,j,jt,rotvec(3),dumint1,jj,jsig
 real*8 dwei,emin,emax,rvec(3),btransform(3,3),dw
 integer ndim(3),nnewsk,iarg,kvec(3),count_lines,kdesym(4,3)
 integer,allocatable :: newmap(:,:)
 character*6 inputkgen
 character*1 dum1
 character*80 klistfilename,klistfullfilename
 logical lattconv,inter
 integer unitstruct,unitkgen,unitinput,nk2
 type(structure) lattice

 iarg=iargc() 
 theta = 1d0
 roundi = 10 !number of digits the mapping is rounded to
 inter=.false.
!  rounder = 100000000000.0 
 lattconv = .false.
 if(iarg.ge.1) then
    do j=1,iarg
       call getarg(j,argdummy)
       if (argdummy(1:1).eq.'-') then
             if ((argdummy(2:2).eq.'h').or.(argdummy(2:2).eq.'H')) then
                  write(*,*)"adaptive refinement of tetrahedral k-mesh"
                  write(*,*)"input: case.klist: starting symmetrized k-points"
                  write(*,*)"       case.klist_full: starting unsymmetrized k-points"
                  write(*,*)"       case.kgen: starting symmetrized tetrahedra"
                  write(*,*)"       case.kgen_full: starting unsymmetrized tetrahedra"
                  write(*,*)"       case.kcontribw: function values for estimator(referring to klist)"
                  write(*,*)"       case.voe: internal list of k-points on tetrahedral edges"
                  write(*,*)"       case.map: internal mapping of klist_full to klist"
                  write(*,*)"output: all input files that are updated and obtain an additional"
                  write(*,*)"        _refined, that is case.klist_refined, case.kgen_refined..."
                  write(*,*)"options: -th=<value> defines the 'harshness' of the refinement"
                  write(*,*)"            value=0 means uniform refinement, value=1 no refinement"
                  write(*,*)"            0<value<1 means adaptive refinement with decreasing harshness"
                  write(*,*)"         -input=init<value> triggers initial refinement to " 
                  write(*,*)"            provide the overall starting mesh for the algorithm." 
                  write(*,*)"            value defines the number of initial refinements"
                  write(*,*)"            and should be between 2-4"  
                  write(*,*)"         -l leads to klist treatment as in H,R lattices"  
                  write(*,*)"         -inter leads to enhanced refinement for the tetrahedra"
                  write(*,*)"            responsible for transitions in the higher energy region"    
                  stop
             elseif ((argdummy(2:4).eq.'th=')) then     
                  read(argdummy(5:8),"(F4.2)")theta
             elseif ((argdummy(2:7).eq.'input=')) then     
                  read(argdummy(8:11),"(A6)")inputkgen
                  if (inputkgen.eq."init") then
                     read(argdummy(12:13),"(I2)")initsteps
                  endif
             elseif ((argdummy(2:2).eq.'l')) then     
                 lattconv = .true.
             elseif ((argdummy(2:6).eq.'inter')) then     
                 inter = .true.
             else
                  write(*,*) 'Usage: refine_tmesh [-th=<val>,-input=init<val>]case'
                   stop
              endif
       else
               read(argdummy,*)seedname
       endif
!      
      
   enddo
 endif


  unitstruct = 2
  unitkgen = 3
  unitinput = 4
  open(unit=unitinput,file=clearspace(seedname)//'.woptin',status='old')
  read(unitinput,*)
  read(unitinput,*)rdum,dw
  close(unitinput)

  open(unit=unitstruct,file=clearspace(seedname)//'.struct',status='old')
  open(unit=unitkgen,file=clearspace(seedname)//'.outputkgen_orig',status='old')
  call countatoms(unitstruct,lattice)
  call readin_lattice(unitstruct,unitkgen,lattice)
  close(unitstruct)
   !round transform matrix
   do jd1=1,3
     do jd2=1,3
        lattice%transform_matrix(jd1,jd2) = round_real(lattice%transform_matrix(jd1,jd2),4)
     enddo
  enddo
  write(*,*)"  transformation matrix between wien2k and internal representation"
  write(*,"(3F12.5)")lattice%transform_matrix(1,:)
  write(*,"(3F12.5)")lattice%transform_matrix(2,:)
  write(*,"(3F12.5)")lattice%transform_matrix(3,:)

  nvoe = 0

  if (inputkgen.eq."init") then
    !get initial mesh
    nt = 48*8**(initsteps-1) 
    allocate(tvariance(nt),tclass(nt),newtclass(8*nt))
    nkfull = (3*2**initsteps)**3
    allocate(tetra(nt,10),wtetra(nt),newk(nkfull,3),newtetra(48*8**(initsteps-1),10),&
         newwtetra(48*8**(initsteps-1)),trefine(48*8**(initsteps-1)))
    allocate(VOE((3*2**initsteps)**3,300,2),VOEidx(6*48*8**(initsteps),3))
    VOE = 0
    VOEidx = 0
    call get_initialmesh(initsteps,tcount,newtetra,newwtetra,newtclass,nkfull,newk,nvoe,VOE,VOEidx)
    ndim = 2**(initsteps+1)
    nnewk = 0
    kdim = 2**initsteps
    nk = 0
  else
      !standard run: load k-mesh
      open(17,file=trim(seedname)//'.kgen_full')
      read(17,*)dum1,nt,ndim
      allocate(tetra(nt,10),wtetra(nt),newk(48*nt,3),newtetra(8*nt,10),newwtetra(8*nt),trefine(nt))
      allocate(tvariance(nt),tclass(nt),newtclass(8*nt))
      do jt=1,nt
         read(17,*)tetra(jt,:),wtetra(jt),tclass(jt)
      enddo
      close(17)
      klistfilename = trim(seedname)//'.klist'
      klistfullfilename = trim(seedname)//'.klist_full'
      open(17,file=trim(klistfilename))
      nk = line_count(17) - 2 !take care
      allocate(k(nk,3),nabk(nk))
      do j=1,nk
         if(j.eq.1) then 
            if (inputkgen.eq."wien2k") then
                  read(17,1523) IK,(k(j,i),i=1,3),kdiv, &
                           dwei,emin,emax,nkp,ndim
            else
                  read(17,1523) IK,(k(j,i),i=1,3),kdiv, &
                           dwei,emin,emax,nkp,kdim
            endif
         else
         read(17,1520) IK, (k(j,i),i=1,3),kdiv,dwei
         endif 
      enddo
      close(17)
      open(17,file=trim(klistfullfilename))
      nkfull = line_count(17) - 2 !take care
      allocate(kfull(nkfull,3),kfulltmp(nkfull,3),map(nkfull,2),kfull_int(nkfull,3))
      do j=1,nkfull
         if(j.eq.1) then 
         read(17,1523) IK,(kfull(j,i),i=1,3),kdiv, &
                        dwei,emin,emax,nkp,kdim
         else
         read(17,1520) IK, (kfull(j,i),i=1,3),kdiv,dwei
         endif 
      enddo
      close(17)
      1523 FORMAT(I10,4I10,3f5.1,4x,i6,10x,3i3,1x) 
      1520 FORMAT(I10,4I10,f5.1)     

      !map first to internal coordinates(real numerators) and than
      !back to integer values
      do j=1,nkfull
!          write(*,*)"kfull before map",j,kfull(j,:)
         do jd1=1,3
            rvec(jd1) = real(kfull(j,jd1))/real(ndim(jd1))
         enddo
!          write(*,*)"kfull real",j,rvec
         if ((lattice%specifier.ne."R").and.(lattice%specifier.ne."H").and.(lattice%specifier.ne."C").and.(.not.lattconv)) then
            rvec = matmul(lattice%transform_matrix,rvec)
         endif
!          write(*,*)"kfull transformed",j,rvec
!          write(*,*)"kfull w2k internal",j,rvec
         do jd1=1,3
            kfull(j,jd1) = int(rvec(jd1)*ndim(jd1))
         enddo
!          write(*,*)"kfull after map",j,kfull(j,:)
      enddo
      
      !load mesh information
      open(17,file=trim(seedname)//".voe")
      read(17,*)dum1,nvoe
      allocate(VOE(48*nt,300,2),VOEidx(nvoe + 48*nt,3))
      VOE = 0
      VOEidx = 0
      do jv=1,nvoe
        read(17,*)VOEidx(jv,:)
        nev = 1
        do while (VOE(VOEidx(jv,1),nev,1).ne.0)
          nev = nev + 1
        enddo
        VOE(VOEidx(jv,1),nev,1) = VOEidx(jv,2)
        VOE(VOEidx(jv,1),nev,2) = VOEidx(jv,3)
     enddo
     do jv=1,nkfull
       nev = 0
       do while (VOE(jv,nev+2,1).ne.0)
          nev = nev + 1
       enddo
       VOE(jv,1,1) = nev
     enddo
     close(17)
  endif  


 if ((trim(inputkgen).ne."debu").and.(trim(inputkgen).ne."init")) then
    write(*,*)"standard mode"
    open(102,file=trim(seedname)//'.kcontribw')
    read(102,*)dum1,dumint1,nw,nsig
    if (dumint1.ne.nk) then
       write(*,*)"Error in refine_tmesh: number of k-points inconsistent, check .kcontribw and .klist"
       stop
    endif
    allocate(kcontribw(nk,0:nw,nsig),tmp(nsig))
    do jk=1,nk
         read(102,3000)dumint1,(tmp(jj),jj=1,nsig)
         do jsig=1,nsig
             kcontribw(jk,0,jsig) = tmp(jsig)  
         enddo
         do jw=1,nw
           read(102,3001)tmp(1:nsig)
           do jsig=1,nsig
             kcontribw(jk,jw,jsig) = tmp(jsig)  
           enddo
         enddo
    enddo
    close(102)
3000 FORMAT(I8,6E20.12)
3001 FORMAT(6E20.12)
 elseif (trim(inputkgen).eq."debu") then
   write(*,*)"debug mode"
   nsig = 1 
   nw = 0
   allocate(kcontribw(nk,0:nw,nsig),tmp(nsig))
   do jk=1,nk
      kcontribw(jk,0,1) = dexp(dble(k(jk,1)))
    enddo
   kcontribw(1,0,1) = 0d0
    
!    kcontrib = 1d0
!    kcontrib(10) = 2d0
!    kcontrib(20) = 4d0
 else
     nw = 0
     nsig = 1
     allocate(kcontribw(nk,0:nw,nsig),tmp(nsig))
     kcontribw = 1d0
 endif
  sumw = 0d0
  do jt=1,nt
      sumw = sumw + wtetra(jt)
  enddo

  if ((trim(inputkgen).ne."init")) then
      write(*,*)">>> refinement of tetrahedral mesh"
      write(*,*)"  theta=",theta
      write(*,*)"  number of initial tetrahedra:",nt
      write(*,*)"  number of initial full mesh k-points:",nkfull
      write(*,*)"  number of initial symmetrized k-points:",nk
      write(*,*)"  sum of tetrahedron weights before refinement:",sumw
     
      !get symmetry mesh information  
      kfulltmp = kfull
      !call symmetrize_kmesh(0,nkfull,kfulltmp,lattice%nop,lattice%sym,ndim,nk2,map)
      open(17,file=trim(seedname)//'.map')
      do jk=1,nkfull
          read(17,*)map(jk,:)
      enddo
      close(17)
  endif

  
 
    if (inputkgen.ne."init") then
      !compute error estimators and mark elements
      call nullestimator(nt,nkfull,nk,tetra,tclass,wtetra,map,nw,nsig,&
                                   kcontribw,theta,trefine,inter,dw)
      flush(6) 
      newk(1:nkfull,:) = kfull !remember old k points
      !refine marked elements
      call refine_tetra(nt,tetra,tclass,wtetra,nkfull,48*nt,newk,&
                           nvoe,VOE,VOEidx,trefine,8*nt,&
                            newtetra,newwtetra,newtclass,nnewk,tcount) 
      ndim = ndim*2 !doubled divisor for new mesh
  endif

 !write-out
 open(17,file=trim(seedname)//".voe_refined")
 write(17,*)"# ",nvoe
 do jv=1,nvoe
     write(17,*)VOEidx(jv,:)
!      VOE(VOEidx(jv,1),VOEidx(jv,2)) = VOEidx(jv,3)
 enddo
 close(17)

 !write-out full k mesh and tetrahedra
 if ((lattice%specifier.ne."R").and.(lattice%specifier.ne."H").and.(lattice%specifier.ne."C").and.(.not.lattconv)) then
     call inverse3x3(lattice%transform_matrix,btransform)
     
       do jd1=1,3
           do jd2=1,3
              btransform(jd1,jd2)=round_real(btransform(jd1,jd2),3)
           enddo
       enddo
 else 
    btransform = 0d0
    btransform(1,1) = 1d0
    btransform(2,2) = 1d0
    btransform(3,3) = 1d0
 endif
 write(*,*)"  transformation matrix between internal and wien2k representation"
 write(*,"(3F12.5)")btransform(1,:)
 write(*,"(3F12.5)")btransform(2,:)
 write(*,"(3F12.5)")btransform(3,:)

 !compute the shapes of the tetrahedra
 call compute_shapeparameters(nkfull+nnewk,tcount,&
                      newk(1:nkfull+nnewk,:),newtetra(1:tcount,:),newtclass(1:tcount),ndim)


 open(102,file=trim(seedname)//'.klist_full_refined')
 do jk=1,nkfull+nnewk
   !transform to conventional reciprocal lattice vectors
   rvec = matmul(btransform,dble(newk(jk,:)))
!    write(*,*)jk,rvec
   do jd1=1,3
      !round the result
      rvec(jd1) = round_real(rvec(jd1),2)
   enddo
   if (jk.eq.1) then
       write(102,1533) jk,(int(rvec(ir)),ir=1,3),ndim(1), &
                        rdum,-7.,1.5,nkfull+nnewk,kdim
   else
     write(102,1530) jk, (int(rvec(ir)),ir=1,3),ndim(1),rdum
   endif
 enddo
 write(102,1521)
 close(102)

 open(102,file=trim(seedname)//'.kgen_full_refined')
 write(102,*)"#",tcount,ndim 
 do jt=1,tcount
    write(102,"(10I10,E20.12,1I10)")(newtetra(jt,ir),ir=1,10),newwtetra(jt),newtclass(jt)
 enddo
 close(102)

!  write(*,*)"  number of unsymmetrized k-points ",nkfull
 flush(6)

 !symmetrize kmesh
 allocate(newmap(nkfull+nnewk,2))
 newmap = 0
 if (inputkgen.ne."init") then
    newmap(1:nkfull,:) = map !remember the old map
    k = 2*k !going to new integer mesh with doubled divisor
 endif
 write(*,*)"  symmetrize k-point mesh "
 call symmetrize_kmesh(nk,nkfull+nnewk,newk(1:nkfull+nnewk,:),k,lattice%nop,lattice%ksym,ndim,nnewsk,newmap)

 write(*,*)"  number of full mesh k-points after refinement: ",nkfull+nnewk
 write(*,*)"  number of symmetrized k-points after refinement: ",nnewsk
 write(*,*)"  symmetrize tetra mesh "
 flush(6)

 !symmetrize tetrahedra
 call symmetrize_tetra(nkfull+nnewk,newmap,tcount,newtetra(1:tcount,:),newwtetra(1:tcount),nnewt)
 write(*,*)"  number of full mesh tetrahedra after refinement: ",tcount
 write(*,*)"  number of symmetrized tetrahedra after refinement: ",nnewt
 sumw = 0d0
 do jt=1,nnewt
   sumw = sumw + newwtetra(jt)
 enddo
 write(*,*)"  sum of tetrahedron weights after refinement",sumw
  
 !write-out symmetrized mesh
 open(102,file=trim(seedname)//'.klist_refined')
 do jk=1,nnewsk
   !transform to conventional reciprocal vectors,
   !but only the new ones, since the old ones have already been transformed
   if (jk.gt.nk) then
      rvec = matmul(btransform,dble(newk(jk,:)))
      do jd1=1,3
          rvec(jd1) = round_real(rvec(jd1),2)
      enddo
   else
      rvec = dble(k(jk,:))
   endif

   if (jk.eq.1) then
       write(102,1533) jk,(int(rvec(ir)),ir=1,3),ndim(1), &
                        rdum,-7.,1.5,nnewsk,kdim
   else
     write(102,1530) jk, (int(rvec(ir)),ir=1,3),ndim(1),rdum
   endif
 enddo
 write(102,1521)
 close(102)

 !write-out map
 open(102,file=trim(seedname)//'.map_refined')
 do jk=1,nkfull+nnewk
    write(102,*)newmap(jk,:)
 enddo
 close(102)

 open(102,file=trim(seedname)//'.klist_new')
 do jk=nk+1,nnewsk
   !transform to conventional reciprocal lattice vectors
   rvec = matmul(btransform,dble(newk(jk,:)))
   do jd1=1,3
        rvec(jd1) = round_real(rvec(jd1),2)
   enddo
   if (jk.eq.nk+1) then
       write(102,1533) jk,(int(rvec(ir)),ir=1,3),2*kdiv, &
                        rdum,-7.,1.5,nnewsk,kdim
   else
     write(102,1530) jk, (int(rvec(ir)),ir=1,3),2*kdiv,rdum
   endif
 enddo
 write(102,1521)
 close(102)

 open(102,file=trim(seedname)//'.kgen_refined')
 write(102,*)"#",nnewt,ndim
 do jt=1,nnewt
    write(102,"(10I10,E20.12)")(newtetra(jt,ir),ir=1,10),newwtetra(jt) 
 enddo
 close(102)
 
 1530 FORMAT(I10,4I10,f5.1)                                              
 1521 format('END',/)
 1533 FORMAT(I10,4I10,3f5.1,4x,i6,' k, div: (',3i3,')')    
 
 CONTAINS
 FUNCTION kmean(kp1,kp2,ndim) result(outm)
    
  implicit none

  integer kp1(3),kp2(3),ndim(3),outm(3),jd

  outm = 1
  do jd=1,3
     if (kp1(jd).gt.ndim(jd)) then
        kp1(jd) = 2*ndim(jd)-kp1(jd)
     endif
     if (kp2(jd).gt.ndim(jd)) then
        kp2(jd) = 2*ndim(jd)-kp2(jd)
     endif
     outm(jd) = (kp1(jd)+kp2(jd))/2
  enddo
  END FUNCTION
  
  

 END PROGRAM

  SUBROUTINE desym(nk,k,tetra,ndim,outk) 

  implicit none

  integer nk,k(nk,3),tetra(4),ndim(3),outk(4,3)
  integer jv,jd,idx1(4),idx2(4),jt
  
  outk = 0
!   tmpk = 0
  do jt=1,4
     outk(jt,:) = k(tetra(jt),:)
!     write(*,*)"my",tetra(jt),k(tetra(jt),:)
  enddo
  do jd=1,3
     idx1 = 0
     idx2 = 0
     do jt=1,4
        if (outk(jt,jd).gt.ndim(jd)/2) then
           idx1(jt) = 1
        elseif (outk(jt,jd).eq.0) then
           idx2(jt) = 1
        endif
     enddo
    if ((any(idx1.eq.1)).and.(any(idx2.eq.1))) then
       do jt=1,4
          if (idx2(jt).eq.1) then
             outk(jt,jd) = ndim(jd)
          endif
       enddo
    endif
  enddo

  END SUBROUTINE


!     *****************************************************************
      FUNCTION find_kpoint2(nk,k,kvec) result(outp)

      implicit none

      integer outp,nk,k(max(nk,1),3),kvec(3),j,kmax(3)

      if (nk.eq.0) then
        outp = -1
        return
      endif
      outp = 0
      do j=1,nk
           if ((k(j,1).eq.kvec(1)).and.(k(j,2).eq.kvec(2)).and.(k(j,3).eq.kvec(3))) then
              outp = j
              return
           endif
      enddo

      END FUNCTION


      


  SUBROUTINE symmetrize_kmesh2(nk,k,nop,symop,Ndim,nsk,map)
  !on input k: unsymmetrized klist
  !on output k: symmetrized klist with rest zeros
  
  implicit none

  integer :: nk,nop,jk,js,jd1,jd2,jsk
  integer :: k(nk,3),tmpk(nk,3),nsk
  real*8 :: symop(nop,3,3)
  integer :: rotvec(3),locvec(3),ndim(3)
  integer :: map(nk)
!      write(*,*)"ndim",ndim,nk
  tmpk = 0
  nsk = 0
!   write(*,*)"M",map(1)
!   
!       do jk=1,nk
!           write(*,*)"lock",nk,jk,k(jk,:)
!       enddo
!     stop
  
  do jk=1,nk
      write(*,*)"lock",ndim(1),k(jk,:)
    locvec = k(jk,:)
    do js=1,nop
         rotvec = 0
         do jd1=1,3
            do jd2=1,3
               rotvec(jd1) = rotvec(jd1) + int(symop(js,jd1,jd2))*k(jk,jd2)
            enddo
         enddo
         if ((rotvec(1).ge.0).and.(rotvec(2).ge.0).and.(rotvec(3).ge.0)) then
              if (rotvec(1).gt.Ndim(1)/2) then
                  rotvec(1) = Ndim(1) - rotvec(1) 
             endif
             if (rotvec(2).gt.Ndim(2)/2) then
                 rotvec(2) = Ndim(2) - rotvec(2)
             endif
             if (rotvec(3).gt.Ndim(3)/2) then
                rotvec(3) = Ndim(3) - rotvec(3) 
             endif
             if ((rotvec(1).lt.locvec(1)).or.&
                  ((rotvec(1).eq.locvec(1)).and.(rotvec(2).le.locvec(2)))) then
!                 write(*,*)"replace",js,rotvec
               locvec = rotvec
!           write(*,*)"rotvec",js,rotvec
            endif
         endif
     enddo
     k(jk,:) = locvec
!      write(*,*)"finvec",locvec
  enddo
  write(*,*)symop(1,1,:)
  write(*,*)symop(1,2,:)
  write(*,*)symop(1,3,:)

  map = 0
  do jk=1,nk
    do jsk=1,nsk
        if  ((k(jk,1).eq.tmpk(jsk,1)).and.(k(jk,2).eq.tmpk(jsk,2))&
           .and.(k(jk,3).eq.tmpk(jsk,3))) then
          map(jk) = jsk
        endif
    enddo
    if (map(jk).eq.0) then
       nsk = nsk + 1
       tmpk(nsk,:) = k(jk,:)
       map(jk) = nsk
!          write(*,*)"newk",k(jk,:)
     else
!          write(*,*)"mapped",jk,map(jk)
     endif
  enddo

  k = 0
  do jsk=1,nsk
    k(jsk,:) =  tmpk(jsk,:)
!      write(*,*)"ksym",k(jsk,:)
  enddo
  

  END SUBROUTINE

  SUBROUTINE symmetrize_kmesh(nskold,nk,k,skold,nop,symop,Ndim,nsk,map)
  !on input k: unsymmetrized klist
  !on output k: symmetrized klist with rest zeros
  !the plan is to use the symmetry operations on the full kmesh
  !than order the kpoints and then take the first k point
  ! as representative for the symmetrized list
  ! It is thus necessary to remember the original index during the ordering
  implicit none

  integer :: nk,nop,jk,js,jd1,jd2,jsk
  integer :: k(nk,3),tmpk(nk*nop,5),tmpk2(nk*nop,5),nsk
  integer :: symop(nop,3,3)
  integer :: rotvec(3),locvec(3),ndim(3),skold(nskold,3)
  integer :: map(nk,2),mapsmap(nk),ichang,compvec,it(4),counter,jc,jd,jk2
  integer :: idx(nk*nop),flux(ndim(1)+2),flux2(ndim(1)+2),fidx,fidx2,tmp(nk*nop),nf,nf2
  integer :: tmpk3(nk*nop,5),lock,nskold,nkoldcount,nknewcount,nkoldtmp,mapidx(nk*nop),mapcount
  logical :: isold(nk)
!      write(*,*)"ndim",ndim,nk,nskold,nop

  write(*,*)"  number of symmetry operations:",nop
!   do js=1,nop
!      do jd1=1,3
!         write(*,*)js,symop(js,jd1,:)
!      enddo
!   enddo
  tmpk = 0
  nsk = 0
  isold = .false.

!     do jk=1,nk
!         write(*,*)"klist before sym",nk,jk,k(jk,:)
!     enddo
!     stop
!    do js=1,nop
!       symop(js,1,:) = (/1d0,0d0,0d0/)
!       symop(js,2,:) = (/0d0,1d0,0d0/)
!       symop(js,3,:) = (/0d0,0d0,1d0/)
!    enddo

  counter =1
  tmpk = 0
  !act all symmetry operations on all unsymmetrized k-points
  do js=1,nop
      do jk=1,nk
         rotvec = 0   
         do jd1=1,3
            do jd2=1,3
               rotvec(jd1) = rotvec(jd1) + 2*int(symop(js,jd1,jd2))*k(jk,jd2)
            enddo
            rotvec(jd1) = mod(rotvec(jd1),2*Ndim(jd1))
            rotvec(jd1) = rotvec(jd1) + (1-isign(1,rotvec(jd1)))*Ndim(jd1)
         enddo
            
         tmpk(counter,1:3) = rotvec(1:3)/2
         tmpk(counter,4) = jk
         tmpk(counter,5) = js
         counter = counter + 1
      enddo
  enddo

  !order according to first dimensional entry
  call indexx(tmpk,idx,nop*nk)
    tmpk2 = 0
  do jk=1,nop*nk
    tmpk2(jk,:) = tmpk(idx(jk),:)
!       write(*,*)"s1",jk,tmpk2(jk,:)
  enddo

  !now go for the other entries
  !the klist now has intervals of equal 1st entries
  !which have to be ordered according to the second entry
  flux(1) = 0
  fidx = 1
  do jk=2,nop*nk
     !find the number of intervals and their boundaries
     if (tmpk2(jk-1,1).ne.tmpk2(jk,1)) then
        fidx = fidx + 1
        flux(fidx) = jk - 1
!         write(*,*)flux(fidx)
     endif
  enddo
  flux(fidx+1) = nop*nk
  tmpk = 0
  tmpk3 = 0
  do jk=2,fidx+1
     nf = flux(jk) - flux(jk-1)
     tmp(1:nf) = tmpk2(flux(jk-1)+1:flux(jk),2)
     call indexx(tmp,idx,nf)
     tmpk(flux(jk-1)+1:flux(jk),:) = tmpk2(flux(jk-1)+idx(1:nf),:)

     !now the same game for the third entry
     fidx2 = 1
     flux2 = 0
     flux2(1) = flux(jk-1) 
     do jk2=flux(jk-1)+2,flux(jk)
        if (tmpk(jk2-1,2).ne.tmpk(jk2,2)) then
          fidx2 = fidx2 + 1
          flux2(fidx2) = jk2 - 1
        endif
     enddo
     flux2(fidx2+1) = flux(jk)
     do jk2=2,fidx2 + 1
        nf2 = flux2(jk2) - flux2(jk2-1)
        tmp(1:nf2) = tmpk(flux2(jk2-1)+1:flux2(jk2),3)
        call indexx(tmp,idx,nf2)
        tmpk3(flux2(jk2-1)+1:flux2(jk2),:) = tmpk(flux2(jk2-1)+idx(1:nf2),:)
!          do jd1=(flux2(jk2-1)+1),flux2(jk2)
!              write(*,*)"s",tmpk3(jd1,:)
!          enddo
     enddo
  enddo

  nsk = nskold

  k(1:nsk,:) = skold
  tmpk = 0
  tmpk(1,:) = tmpk3(1,:)
  lock = 1

  mapcount = 1
  mapidx = 0 
  mapidx(1) = map(tmpk3(1,4),1)!map(tmpk3(1,4),1)

  !now take a first look at the ordered list
  !and check whether the old symmetrized klist 
  !already covers a k point
  do jk = 2,nop*nk
      if ((tmpk3(jk,1).ne.tmpk3(jk-1,1)).or.&
          (tmpk3(jk,2).ne.tmpk3(jk-1,2)).or.&
          (tmpk3(jk,3).ne.tmpk3(jk-1,3))) then
          mapcount = mapcount + 1
      endif
      if (mapidx(mapcount).eq.0) then
           mapidx(mapcount) = map(tmpk3(jk,4),1)
      endif
  enddo


  !here comes the real mapping
  if (mapidx(1).eq.0) then
      nsk = nsk + 1
      map(tmpk3(1,4),1) = nsk
      mapidx(1) = nsk
  else
      map(tmpk3(1,4),1) = mapidx(1)
  endif
  map(tmpk3(1,4),2) = tmpk3(1,5)
  mapcount = 1
  do jk = 2,nop*nk
!        write(*,*)"kpoint",jk,map(tmpk3(jk,4),1),tmpk3(jk,:)
      if ((tmpk3(jk,1).ne.tmpk3(jk-1,1)).or.&
          (tmpk3(jk,2).ne.tmpk3(jk-1,2)).or.&
          (tmpk3(jk,3).ne.tmpk3(jk-1,3))) then
          mapcount = mapcount + 1
!              write(*,*)"new interval in klist",jk,mapcount
      endif
      if (map(tmpk3(jk,4),1).eq.0) then
        !k-point has no mapping
	if (mapidx(mapcount).eq.0) then
               !no other k-point in this symmetry group 
               !has a mapping -> it has to be a new one
		nsk = nsk + 1
	        k(nsk,:) = tmpk3(jk,1:3)
		mapidx(mapcount) = nsk
		map(tmpk3(jk,4),1) = nsk
                map(tmpk3(jk,4),2) = tmpk3(jk,5)
	else
                !another k-point in this symmetry group 
                !has a mapping -> map to a kpoint from the old symmetrized klist 
		map(tmpk3(jk,4),1) = mapidx(mapcount)
                map(tmpk3(jk,4),2) = tmpk3(jk,5)
	endif
      endif
  enddo
   
   do jk=1,nk
      if (map(jk,1).eq.0) then
          write(*,*)"Error in symmetrize k-mesh",jk
          stop
      endif
   enddo

  END SUBROUTINE

  SUBROUTINE symmetrize_tetra(maxk,map,nt,tetra,wtetra,countsymtetra)

  implicit none

  integer maxk,map(maxk),nt,tetra(nt,10),compvec,jc
  real*8 wtetra(nt),tmpw(nt),r1
  integer jt,jd,jd1,jd2,val1,val2,idx(nt),tmpt(nt,10)
  integer :: i1,i2,i3,i4,ichang,loctetra(10),it(10)
  integer :: countsymtetra

  !first map the tetrahedra vertices to symmetrized kmesh
  do jt=1,nt
     do jd=1,10       
       tetra(jt,jd) = map(tetra(jt,jd))
     enddo
  enddo

  !now order vertices for each tetrahedron from lowest to highest value
  do jt=1,nt
      do jd1=1,3                                                       
         do jd2=jd1+1,4                                                     
            val1=tetra(jt,jd1)                                                  
            val2=tetra(jt,jd2)                                                  
            tetra(jt,jd1)= min(val1,val2)                                     
            tetra(jt,jd2)= max(val1,val2)
         enddo
      enddo
  enddo
  do jt=1,nt
      do jd1=5,9                                                       
         do jd2=jd1+1,10                                                     
            val1=tetra(jt,jd1)                                                  
            val2=tetra(jt,jd2)                                                  
            tetra(jt,jd1)= min(val1,val2)                                     
            tetra(jt,jd2)= max(val1,val2)
         enddo
      enddo
  enddo

  !reorder tetra list according to the first vertex
  !heapsort from wien2k
  call indexx(tetra,idx,nt)
  do jt=1,nt
    tmpt(jt,:) = tetra(idx(jt),:)
    tmpw(jt) = wtetra(idx(jt))
  enddo


do jd=2,10
   ichang = 1
   do while (ichang.eq.1)
         ichang = 0
         do jt=2,nt
            compvec = 0
            do jc=1,jd-1
              if (tmpt(jt,jc).eq.tmpt(jt-1,jc)) then
                 compvec = compvec + 1
              endif 
            enddo
            if ((compvec.eq.jd-1).and.(tmpt(jt,jd).lt.tmpt(jt-1,jd))) then
                do jc=jd,10
                  it(jc)=tmpt(jt-1,jc)   
                  tmpt(jt-1,jc)=tmpt(jt,jc)
                  tmpt(jt,jc)=it(jc)
                enddo
                r1=tmpw(jt-1)
                tmpw(jt-1)=tmpw(jt)
                tmpw(jt)=r1
                ichang=1
            end if
         enddo
   enddo
enddo

   !count the tetrahedra and sum up the weight
   countsymtetra = 1
   jt=1
   loctetra = tmpt(1,:)
   tetra = 0
   wtetra = 0d0
   do while (jt.le.nt)
     if ((tmpt(jt,1).eq.loctetra(1)).and.(tmpt(jt,2).eq.loctetra(2)).and.&
         (tmpt(jt,3).eq.loctetra(3)).and.(tmpt(jt,4).eq.loctetra(4)).and.& 
         (tmpt(jt,4).eq.loctetra(4)).and.(tmpt(jt,5).eq.loctetra(5)).and.&
         (tmpt(jt,6).eq.loctetra(6)).and.(tmpt(jt,7).eq.loctetra(7)).and.&
         (tmpt(jt,8).eq.loctetra(8)).and.(tmpt(jt,9).eq.loctetra(9)).and.&
         (tmpt(jt,10).eq.loctetra(10))) then
         wtetra(countsymtetra) = wtetra(countsymtetra) + tmpw(jt)
     else
         tetra(countsymtetra,:) = loctetra
         countsymtetra = countsymtetra + 1
         loctetra = tmpt(jt,:)
         wtetra(countsymtetra) = tmpw(jt)
                   
     endif
     jt=jt+1
   enddo
   tetra(countsymtetra,:) = loctetra
   
     
  END SUBROUTINE

 FUNCTION test1(k) result(outk)
 
 implicit none

 integer k,outk

!  write(*,*)"b1",k
 if (k.gt.16) then
   outk = 32 - k
 else
    outk = k
 endif
!  write(*,*)"b2",outk

 END FUNCTION
 
 SUBROUTINE nullestimator(nt,nkfull,nk,tetra,tclass,wtetra,map,nw,nsig,kcontribw,theta,trefine,inter,dw)


  implicit none

  integer :: nt,nkfull,nk,jt,jd,jd2,jk,idxfull(nkfull),counterm,nkfree,nkfree2,tclass(nt)
  integer :: tetra(nt,10),map(nkfull),info,jk2,test1,nw,nsig,jsig,jw,reftet(8,4)
  integer :: jref
  real*8 :: kcontribw(nk,0:nw,nsig),estimator(nt),maxest,theta,sumest
  real*8 :: wtetra(nt),est1,est2,est3,contrw,dw!,mmatrix(nkfull,nkfull),!refmat(4,4)!,RHS(nkfull)
  logical :: trefine(nt),inter
!   integer, allocatable ::  idx(:),ipiv(:)
!   real*8, allocatable :: mmatrixred(:,:),RHSred(:)
  
  
  write(*,*)" error estimators for frequencies [w errorestimator]"
  estimator = 0d0
  do jw=0,nw
     contrw = 0d0
     do jsig=1,nsig
         do jt=1,nt
            reftet(1,1:4) = (/tetra(jt,1),tetra(jt,5),tetra(jt,6),tetra(jt,7)/)
            reftet(2,1:4) = (/tetra(jt,2),tetra(jt,5),tetra(jt,8),tetra(jt,9)/)
               reftet(3,1:4) = (/tetra(jt,3),tetra(jt,6),tetra(jt,8),tetra(jt,10)/)
               reftet(4,1:4) = (/tetra(jt,4),tetra(jt,7),tetra(jt,9),tetra(jt,10)/)
               if (tclass(jt).eq.1) then !central octahedron
                     !corresp. to class 1 in Endres paper
                     reftet(5,1:4) = (/tetra(jt,5),tetra(jt,6),tetra(jt,7),tetra(jt,10)/)
                     reftet(6,1:4) = (/tetra(jt,5),tetra(jt,6),tetra(jt,8),tetra(jt,10)/)
                     reftet(7,1:4) = (/tetra(jt,5),tetra(jt,7),tetra(jt,9),tetra(jt,10)/)
                     reftet(8,1:4) = (/tetra(jt,5),tetra(jt,8),tetra(jt,9),tetra(jt,10)/)
               else
                  !corresp. to class 2 in Endres paper
                  reftet(5,1:4) = (/tetra(jt,5),tetra(jt,6),tetra(jt,7),tetra(jt,9)/)
                  reftet(6,1:4) = (/tetra(jt,5),tetra(jt,6),tetra(jt,8),tetra(jt,9)/)
                  reftet(7,1:4) = (/tetra(jt,6),tetra(jt,7),tetra(jt,9),tetra(jt,10)/)
                  reftet(8,1:4) = (/tetra(jt,6),tetra(jt,8),tetra(jt,9),tetra(jt,10)/)
               endif
            est1 = 0d0
            est2 = 0d0
            est3 = 0d0
            do jd=1,4
               est1 = est1 + wtetra(jt)*kcontribw(map(tetra(jt,jd)),jw,jsig)/4d0
               est2 = est2 - wtetra(jt)*kcontribw(map(tetra(jt,jd)),jw,jsig)/20d0
!                 write(*,*)" freq,",jw,jsig,jt,jd,map(tetra(jt,jd)),kcontribw(map(tetra(jt,jd)),jw,jsig),&
!                             wtetra(jt)
            enddo
!              write(*,*)"seq1",est1,est2,est3
            do jd=5,10
               est2 = est2 +&
                  wtetra(jt)*kcontribw(map(tetra(jt,jd)),jw,jsig)/5d0
!                write(*,*)" freq,",jw,jsig,jt,jd,map(tetra(jt,jd)),kcontribw(map(tetra(jt,jd)),jw,jsig)
            enddo
!               write(*,*)"seq2",est1,est2,est3
!             write(*,*)" freqc,",jw,jsig,jt,dsqrt((est1(jt)-est2(jt))**2)/dble(nsig)/dble(nw)
            do jref=1,8
               do jd=1,4
                   est3 = est3 + wtetra(jt)*kcontribw(map(reftet(jref,jd)),jw,jsig)/32d0
               enddo
            enddo
            !write(*,*)"seq3",jw,jsig,est1,est2,est3
            if (inter) then
!                write(*,*)"inter is",inter
               estimator(jt) = estimator(jt) + dabs(est1-est3)/dble(nsig)*jw*dw!**2/dabs(est2-est3)/dble(nsig)
               contrw = contrw + dabs(est1-est3)/dble(nsig)*jw*dw!**2/dabs(est2-est3)/dble(nsig)
            else
!                write(*,*)"inter is",inter
               estimator(jt) = estimator(jt) + dabs(est1-est3)/dble(nsig)!**2/dabs(est2-est3)/dble(nsig)
               contrw = contrw + dabs(est1-est3)/dble(nsig)!**2/dabs(est2-est3)/dble(nsig)
            endif
         enddo
     enddo
     write(*,"(2F12.5)")jw*dw,contrw
  enddo

  !find maximal entry 
  maxest = 0d0
  do jt=1,nt
     if (estimator(jt).gt.maxest) then
        maxest = estimator(jt)
     endif
  enddo

  write(*,*)"  maximal error estimator: ",maxest
  !mark elements to refine
  trefine = .false.
  do jt=1,nt
     if ((estimator(jt).gt.theta*maxest).or.(theta.eq.0d0)) then
        trefine(jt) = .true.
     endif
  enddo

  sumest = 0d0
  write(*,*)" error estimators for elements [elementidx errorestimator marked weight vertices]"
  do jt=1,nt
    sumest = sumest + estimator(jt)
     write(*,"(I10,F13.8,L3,F13.8,4I10)")jt,estimator(jt),trefine(jt),wtetra(jt),tetra(jt,1:4)
  enddo
  write(*,*)"estimator = ",sumest
   
  END SUBROUTINE
 
 SUBROUTINE refine_tetra(nt,tetra,tclass,wtetra,nkfull,maxk,newk,nvoe,VOE,VOEidx,trefine,maxt, &
                           newtetra,newwtetra,newtclass,nnewk,tcount) 

 implicit none

 integer maxt,nt,maxk
 integer nkfull,tetra(nt,10),tclass(nt),newtetra(maxt,10),newtclass(maxt),tcount
 integer tchange,jt,npidx(12,2),VOE(maxk,300,2),nvoe,VOEidx(nvoe + 48*nt,3)
 integer jd,jk,kidx(6),locv,nev,jv,nnewk,newk(maxk,3)
 integer eidx1,eidx2,eidx11,eidx22
 real*8 :: wtetra(nt),newwtetra(maxt)
 logical :: trefine(maxt)

 
!  write(*,*)"notorious",nt,tetra(nt,:)
!  write(*,*)"notorious2",nt,tclass(nt),wtetra(nt),nkfull,newk(nkfull,:)
!  write(*,*)"serious?",nvoe + 48*nt
!  write(*,*)"notorious3",nvoe,VOEidx(1,3)
!  write(*,*)"notorious4",maxt!,newtetra(maxt,10)
!  write(*,*)"14rin",newk(14,:)
! write(*,*)"nllb3",nvoe,nt,nvoe + 8*nt
  tcount = nt
  tchange = 1
  newtetra(1:nt,:) = tetra
  newwtetra(1:nt) = wtetra
  newtclass(1:nt) = tclass
  nnewk = 0
  !update the original k points
   do jk=1,nkfull
       newk(jk,:) = 2*newk(jk,:)
! !         write(*,*)"k",jk,2*kfull(jk,:)
   enddo 
!   !update ndim
! !   ndim = 2*ndim
  npidx(1,:) = (/ 1, 2 /)
  npidx(2,:) = (/ 1, 3 /)
  npidx(3,:) = (/ 1, 4 /)
  npidx(4,:) = (/ 2, 3 /)
  npidx(5,:) = (/ 2, 4 /)
  npidx(6,:) = (/ 3, 4 /)

!  do jt=1,nkfull
!   write(*,*)"nllb",VOE(jt,1:10,1)
!  enddo
! 
 tchange = 1
 do while (tchange.ne.0)
   tchange = 0
   do jt=1,nt
!       write(*,*)"tetra0",jt,tetra(jt,:)
      if (trefine(jt)) then
         if ((tetra(jt,1).eq.tetra(jt,2)).or.(tetra(jt,2).eq.tetra(jt,3)).or.(tetra(jt,3).eq.tetra(jt,4))) then
            write(*,*)"Error: multiple vertices",jt,tetra(jt,:)
            stop
         else
               !built 8 new tetrahedron and weights
!pw                newtetra(jt,1:4) = (/tetra(jt,1),tetra(jt,5),tetra(jt,6),tetra(jt,7)/)
!pw                newtetra(tcount+1,1:4) = (/tetra(jt,2),tetra(jt,5),tetra(jt,8),tetra(jt,9)/)
!pw                newtetra(tcount+2,1:4) = (/tetra(jt,3),tetra(jt,6),tetra(jt,8),tetra(jt,10)/)
!pw                newtetra(tcount+3,1:4) = (/tetra(jt,4),tetra(jt,7),tetra(jt,9),tetra(jt,10)/)
               newtetra(jt,1:4) = (/tetra(jt,1),tetra(jt,5),tetra(jt,6),tetra(jt,7)/)
               newtetra(tcount+1,1:4) = (/tetra(jt,5),tetra(jt,2),tetra(jt,8),tetra(jt,9)/)
               newtetra(tcount+2,1:4) = (/tetra(jt,6),tetra(jt,8),tetra(jt,3),tetra(jt,10)/)
               newtetra(tcount+3,1:4) = (/tetra(jt,7),tetra(jt,9),tetra(jt,10),tetra(jt,4)/)
               newtclass(jt) = tclass(jt)
               newtclass(tcount+1) = tclass(jt)
               newtclass(tcount+2) = tclass(jt)
               newtclass(tcount+3) = tclass(jt)
               if (tclass(jt).eq.1) then !central octahedron
                     !corresp. to class 1 in Endres paper
!                       newtetra(tcount+4,1:4) = (/tetra(jt,5),tetra(jt,6),tetra(jt,7),tetra(jt,10)/)
!                       newtetra(tcount+5,1:4) = (/tetra(jt,5),tetra(jt,6),tetra(jt,8),tetra(jt,10)/)
!                       newtetra(tcount+6,1:4) = (/tetra(jt,5),tetra(jt,7),tetra(jt,9),tetra(jt,10)/)
!                       newtetra(tcount+7,1:4) = (/tetra(jt,5),tetra(jt,8),tetra(jt,9),tetra(jt,10)/)
                      newtetra(tcount+4,1:4) = (/tetra(jt,9),tetra(jt,7),tetra(jt,10),tetra(jt,5)/)
                      newtetra(tcount+5,1:4) = (/tetra(jt,10),tetra(jt,6),tetra(jt,5),tetra(jt,7)/)
                      newtetra(tcount+6,1:4) = (/tetra(jt,9),tetra(jt,8),tetra(jt,5),tetra(jt,10)/)
                     newtetra(tcount+7,1:4) = (/tetra(jt,6),tetra(jt,5),tetra(jt,8),tetra(jt,10)/)
                      newtclass(tcount+4) = -1
                      newtclass(tcount+5) = -1
                      newtclass(tcount+6) = 1
                      newtclass(tcount+7) = 1
! 		     newtclass(tcount+4) = 1
!                      newtclass(tcount+5) = -1
!                      newtclass(tcount+6) = 1
!                      newtclass(tcount+7) = -1
               else
                  !corresp. to class 2 in Endres paper
!                     newtetra(tcount+4,1:4) = (/tetra(jt,5),tetra(jt,6),tetra(jt,7),tetra(jt,9)/)
!                     newtetra(tcount+5,1:4) = (/tetra(jt,5),tetra(jt,6),tetra(jt,8),tetra(jt,9)/)
!                     newtetra(tcount+6,1:4) = (/tetra(jt,6),tetra(jt,7),tetra(jt,9),tetra(jt,10)/)
!                     newtetra(tcount+7,1:4) = (/tetra(jt,6),tetra(jt,8),tetra(jt,9),tetra(jt,10)/)
                   newtetra(tcount+4,1:4) = (/tetra(jt,10),tetra(jt,6),tetra(jt,8),tetra(jt,9)/)
                   newtetra(tcount+5,1:4) = (/tetra(jt,10),tetra(jt,9),tetra(jt,7),tetra(jt,6)/)
                   newtetra(tcount+6,1:4) = (/tetra(jt,5),tetra(jt,8),tetra(jt,6),tetra(jt,9)/)
                   newtetra(tcount+7,1:4) = (/tetra(jt,5),tetra(jt,7),tetra(jt,9),tetra(jt,6)/)
                    newtclass(tcount+4) = -1
                    newtclass(tcount+5) = -1
                    newtclass(tcount+6) = 1
                    newtclass(tcount+7) = 1
!                   newtclass(tcount+4) = 1
!                   newtclass(tcount+5) = 1
!                   newtclass(tcount+6) = -1
!                   newtclass(tcount+7) = -1
               endif
               do jd=0,7
                  if (jd.eq.0) then
                     newwtetra(jt) = wtetra(jt)/8d0
                     npidx(1,:) = (/ newtetra(jt,1), newtetra(jt,2)/)
                     npidx(2,:) = (/ newtetra(jt,1), newtetra(jt,3)/)
                     npidx(3,:) = (/ newtetra(jt,1), newtetra(jt,4)/)
                     npidx(4,:) = (/ newtetra(jt,2), newtetra(jt,3)/)
                     npidx(5,:) = (/ newtetra(jt,2), newtetra(jt,4)/)
                     npidx(6,:) = (/ newtetra(jt,3), newtetra(jt,4)/)
                  else
                     newwtetra(tcount+jd) = wtetra(jt)/8d0
                     npidx(1,:) = (/ newtetra(tcount+jd,1), newtetra(tcount+jd,2)/)
                     npidx(2,:) = (/ newtetra(tcount+jd,1), newtetra(tcount+jd,3)/)
                     npidx(3,:) = (/ newtetra(tcount+jd,1), newtetra(tcount+jd,4)/)
                     npidx(4,:) = (/ newtetra(tcount+jd,2), newtetra(tcount+jd,3)/)
                     npidx(5,:) = (/ newtetra(tcount+jd,2), newtetra(tcount+jd,4)/)
                     npidx(6,:) = (/ newtetra(tcount+jd,3), newtetra(tcount+jd,4)/)
                  endif
                     
                  do jk=1,6
                        kidx(jk) = 0
                        eidx1 = min(npidx(jk,1),npidx(jk,2))
                        eidx2 = max(npidx(jk,1),npidx(jk,2))
                        nev = VOE(eidx1,1,1)
                        do jv=1,nev
                           if (VOE(eidx1,jv+1,1).eq.eidx2) then
                              kidx(jk) = VOE(eidx1,jv+1,2)
                              exit
                           endif
                        enddo
                        if (kidx(jk).eq.0) then
                           nnewk = nnewk + 1
                           kidx(jk) = nkfull+nnewk
                           newk(kidx(jk),:) = (newk(npidx(jk,1),:)+newk(npidx(jk,2),:))/2
                           VOE(eidx1,1,1) = nev + 1
                           VOE(eidx1,nev+2,1) = eidx2
!                            write(*,*)"prepare"
!                            write(*,*)"done",eidx1,nev,VOE(eidx1,nev+2,2),kidx(jk)
                           VOE(eidx1,nev+2,2) = kidx(jk)
!                             stop
                           nvoe = nvoe + 1
!                            write(*,*)"getta",jt,jd,jk,nvoe
                           VOEidx(nvoe,:) = (/ min(npidx(jk,1),npidx(jk,2)) ,max(npidx(jk,1),npidx(jk,2)), kidx(jk)/)
!                         else
!                             write(*,*)"old",nev,eidx1,eidx2,kidx(jk)
                        endif
                     if (jd.eq.0) then
                        newtetra(jt,4+jk) = kidx(jk)
                     else
                        newtetra(tcount+jd,4+jk) = kidx(jk)
                     endif
                  enddo
               enddo
               tcount = tcount + 7
         endif
      endif
   enddo
 
    !homogenize mesh
    do jt=1,nt
       if (.not.trefine(jt)) then
         npidx(1,:) = (/ newtetra(jt,1), newtetra(jt,5)/)
         npidx(2,:) = (/ newtetra(jt,1), newtetra(jt,6)/)
         npidx(3,:) = (/ newtetra(jt,1), newtetra(jt,7)/)
         npidx(4,:) = (/ newtetra(jt,2), newtetra(jt,8)/)
         npidx(5,:) = (/ newtetra(jt,2), newtetra(jt,9)/)
         npidx(6,:) = (/ newtetra(jt,3), newtetra(jt,10)/)
         npidx(7,:) = (/ newtetra(jt,2), newtetra(jt,5)/)
         npidx(8,:) = (/ newtetra(jt,3), newtetra(jt,6)/)
         npidx(9,:) = (/ newtetra(jt,4), newtetra(jt,7)/)
         npidx(10,:) = (/ newtetra(jt,3), newtetra(jt,8)/)
         npidx(11,:) = (/ newtetra(jt,4), newtetra(jt,9)/)
         npidx(12,:) = (/ newtetra(jt,4), newtetra(jt,10)/)
         do jk=1,12
            locv = 0
            eidx1 = min(npidx(jk,1),npidx(jk,2))
            eidx2 = max(npidx(jk,1),npidx(jk,2))
!             write(*,*)"stillcounting",eidx1,eidx2
            nev = VOE(eidx1,1,1)
            do jv=1,nev
               if (VOE(eidx1,jv+1,1).eq.eidx2) then
                  locv = VOE(eidx1,jv+1,2)
                  exit
               endif
            enddo
            if (locv.ne.0) then
               eidx11 = min(npidx(jk,1),locv)
               eidx22 = max(npidx(jk,1),locv)
               nev = VOE(eidx11,1,1)
               do jv=1,nev
                  if (VOE(eidx11,jv+1,1).eq.eidx22) then
                     trefine(jt) = .true.
                     tchange = 1
!                     write(*,*)"get down, get down"
                     exit
                  endif
               enddo 
               eidx11 = min(npidx(jk,2),locv)
               eidx22 = max(npidx(jk,2),locv)
               nev = VOE(eidx11,1,1)
               do jv=1,nev
                  if (VOE(eidx11,jv+1,1).eq.eidx22) then
                     trefine(jt) = .true.
                     tchange = 1
!                     write(*,*)"get down, get down"
                     exit
                  endif
               enddo
            endif
          enddo
        else
          trefine(jt) = .false.
        endif
    enddo
 enddo

  
 END SUBROUTINE


 SUBROUTINE get_initialmesh(initsteps,nt,tetra,wtetra,tclass,nkfull,kfull,nvoe,VOE,VOEidx)

 implicit none

 integer initsteps,counter,j1,j2,j3,nkfull,kmax,maxt
 integer kfull((3*2**initsteps)**3,3),Ti(4),krot(4,3)
 integer nt,tetra(48*8**(initsteps-1),10),tclass(48*8**(initsteps-1)),newpoints(6,3),kidx(6)
 integer newtetra(48*8**(initsteps-1),10),newtclass(48*8**(initsteps-1))
 integer VOE((3*2**initsteps)**3,300,2),VOEidx(6*48*8**(initsteps),3)
 real*8 :: S(3,3,48),wtetra(48*8**(initsteps-1)),newwtetra(48*8**(initsteps-1))
 integer eidx1,eidx2,jd,jk,jr,js,jt,nnewk,nvoe,tcount,find_kpoint2,npidx(6,2)
 logical :: trefine(48*8**(initsteps-1))

 kmax=(3*2**initsteps)**3
 maxt = 48*8**(initsteps-1)
 kfull = 0

 counter = 1 
 do j1=1,3
    do j2=1,3
        do j3=1,3
          kfull(counter,:) = (/j1-1,j2-1,j3-1/)  
          counter =  counter + 1
        enddo
    enddo
 enddo
 nkfull = counter - 1

 Ti = (/1,2,5,14/)

 S = 0
 open(102,file=trim('symop'))
 read(102,*)!nop
 do js=1,48
    do jd=1,3
      read(102,"(3F8.5)")S(jd,:,js)
    enddo
    read(102,*)
 enddo
 close(102)

 tetra = 0
 do js=1,48
    krot = 0
    do jt=1,4
       do j1=1,3
          do j2=1,3
             krot(jt,j1) = krot(jt,j1) +int(S(j1,j2,js))*kfull(Ti(jt),j2)-int(S(j1,j2,js))
          enddo
          krot(jt,j1) = krot(jt,j1) + 1
       enddo
   enddo
   do jt=1,4
        do jk=1,nkfull
           if ((krot(jt,1).eq.kfull(jk,1)).and.(krot(jt,2).eq.kfull(jk,2))&
               .and.(krot(jt,3).eq.kfull(jk,3))) then
               tetra(js,jt) = jk
           endif
        enddo
    enddo
    tclass(js) = int(S(1,1,js)*S(2,2,js)*S(3,3,js) + S(1,2,js)*S(2,3,js)*S(3,1,js) &
                 + S(1,3,js)*S(2,1,js)*S(3,2,js) - S(1,1,js)*S(2,3,js)*S(3,2,js) &
                 - S(1,2,js)*S(2,1,js)*S(3,3,js) - S(1,3,js)*S(2,2,js)*S(3,1,js))
 enddo
 
 

 wtetra = 0.125d0/6d0
 nnewk = 0
 VOE = 0
 VOEidx = 0
 nvoe = 0
 npidx(1,:) = (/ 1, 2 /)
 npidx(2,:) = (/ 1, 3 /)
 npidx(3,:) = (/ 1, 4 /)
 npidx(4,:) = (/ 2, 3 /)
 npidx(5,:) = (/ 2, 4 /)
 npidx(6,:) = (/ 3, 4 /)
 nt = 48
 kfull = 2*kfull

 do jt=1,nt
     if ((tetra(jt,1).eq.tetra(jt,2)).or.(tetra(jt,2).eq.tetra(jt,3)).or.(tetra(jt,3).eq.tetra(jt,4))) then
         write(*,*)"Error: multiple vertices",jt,tetra(jt,:)
         stop
     else
       newpoints(1,:) =(kfull(tetra(jt,1),:)+kfull(tetra(jt,2),:))/2
       newpoints(2,:) =(kfull(tetra(jt,1),:)+kfull(tetra(jt,3),:))/2
       newpoints(3,:) =(kfull(tetra(jt,1),:)+kfull(tetra(jt,4),:))/2
       newpoints(4,:) =(kfull(tetra(jt,2),:)+kfull(tetra(jt,3),:))/2
       newpoints(5,:) =(kfull(tetra(jt,2),:)+kfull(tetra(jt,4),:))/2
       newpoints(6,:) =(kfull(tetra(jt,3),:)+kfull(tetra(jt,4),:))/2
     endif
     do jd=1,6
         kidx(jd)= find_kpoint2(nkfull+nnewk,kfull(1:nkfull+nnewk,:),newpoints(jd,:))
         if (kidx(jd).le.0) then
            nnewk = nnewk + 1
            kidx(jd) = nkfull+nnewk
!             write(*,*)"np2",jd,newpoints(jd,:),kidx(jd)
            kfull(kidx(jd),:) = newpoints(jd,:)
            eidx1 = min(tetra(jt,npidx(jd,1)),tetra(jt,npidx(jd,2)))
            eidx2 = max(tetra(jt,npidx(jd,1)),tetra(jt,npidx(jd,2)))
            VOE(eidx1,1,1) = VOE(eidx1,1,1) + 1
            VOE(eidx1,VOE(eidx1,1,1)+2,1) = eidx2
            VOE(eidx1,VOE(eidx1,1,1)+2,2) = kidx(jd)
            nvoe = nvoe + 1
            VOEidx(nvoe,:) = (/eidx1  ,eidx2 , kidx(jd) /)
         endif
         tetra(jt,4+jd) = kidx(jd)
      enddo
 enddo

 

 nkfull = nkfull + nnewk
 nnewk = 0
      
 do jr=1,initsteps-1 
    trefine = .true.
    call refine_tetra(nt,tetra(1:nt,:),tclass(1:nt),wtetra(1:nt),nkfull,kmax,kfull,&
                           nvoe,VOE,VOEidx(1:nvoe + 48*nt,:),trefine,maxt,&
                           newtetra,newwtetra,newtclass,nnewk,tcount) 
    nkfull = nkfull + nnewk
    nt = tcount
    tetra = newtetra
    wtetra = newwtetra
    tclass = newtclass 
 enddo
 
 write(*,*)">>> initial tetrahedral mesh"
 write(*,*)"  number of tetrahedra:",nt
 write(*,*)"  number of full mesh k points:",nkfull
!  do jk=1,nkfull
!     write(*,*)"kfull",jk,kfull(jk,:)
!  enddo
!  do jt=1,nt
!     write(*,*)"tetra",jt,tetra(jt,:)
!  enddo

 END SUBROUTINE

 SUBROUTINE compute_shapeparameters(nk,nt,k,tetra,tclass,ndim)

 implicit none

 integer :: nk,nt,jt,jd1,jd2,jk,ndim(3),counter
 integer :: tetra(nt,10),k(nk,3),tclass(nt)
 real*8 :: vol,crad,qual(nt),a(3),b(3),c(3),quals(nt)
 real*8 :: crossab(3),crossca(3),crossbc(3),tmp(3)

 !do jk=1,nk
 !  write(767,*)k(jk,:)
 !enddo
 !stop
write(767,*)k(1,:)
write(767,*)k(2,:)
write(767,*)k(3,:)
 do jt=1,nt
   
   a = dble(k(tetra(jt,2),:) - k(tetra(jt,1),:))/ndim(1)
   b = dble(k(tetra(jt,3),:) - k(tetra(jt,1),:))/ndim(1)
   c = dble(k(tetra(jt,4),:) - k(tetra(jt,1),:))/ndim(1)
   
   crossab(1) = a(2)*b(3) - a(3)*b(2)
   crossab(2) = a(3)*b(1) - a(1)*b(3)
   crossab(3) = a(1)*b(2) - a(2)*b(1)
   crossca(1) = c(2)*a(3) - c(3)*a(2)
   crossca(2) = c(3)*a(1) - c(1)*a(3)
   crossca(3) = c(1)*a(2) - c(2)*a(1)
   crossbc(1) = b(2)*c(3) - b(3)*c(2)
   crossbc(2) = b(3)*c(1) - b(1)*c(3)
   crossbc(3) = b(1)*c(2) - b(2)*c(1)
   vol = dabs(dot_product(a,crossbc))/6d0
   tmp = dot_product(a,a)*crossbc+dot_product(b,b)*crossca &
        +dot_product(c,c)*crossab
   crad = dsqrt(tmp(1)**2+tmp(2)**2+tmp(3)**2)/12d0/vol
   qual(jt) = crad**3/vol
   write(767,'(A16,I7,3E12.4,I3)')"[idx vol R qual cl]",jt,vol,crad,qual(jt),tclass(jt)
 enddo

 quals = 0d0
 counter = 0
 do jt=1,nt
    if (.not.any(quals(1:counter).eq.qual(jt))) then
       counter = counter + 1
       quals(counter) = qual(jt)
    endif
 enddo

 write(*,*)"  there are ",counter," different shapes"
 write(*,*)"  qualities ",quals(1:counter)

 END SUBROUTINE
