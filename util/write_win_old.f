PROGRAM write_win
 !program prepares input case.win for wannier90 
 
 use util 

 implicit none
 character*50 seedname
 !parameter (seedname='SrVO3')
 integer line_count,nkpoints,j1,j2,tmp
 integer divisor,mp_grid(3)!k-mesh specifiers
 integer iarg,argcount !command line input argument counter
 integer vec(3)
 integer jelement,jmult,jwrite !loop variables
 !integer mult !multiplicity
 integer atomcount !counter for unit cell atoms
 integer nlines   !# lines in a file
 integer nkpoints_band !#kpoints for bandstructure plots
 integer nemin,nemax,tmpint !read in dummys
 integer counter,j
 integer nkpoints_firstpath  !number of kpoints along the first part of the kpath
 real*8, allocatable :: kpoints(:,:),kpoints_band(:,:) !k-mesh points
 real*8 :: rvec(3)
 character*3, allocatable :: kpoints_band_spec(:)
 real*8 :: tmp2(3)
 character*60 dummy
 character*10 dum2
 character*4 dum3!,unit
 character*70 argdummy,startmessage

 logical write_kpath,second_kpath_point

 type(structure) lattice
 type(optiontype) options
 
 character*8 w2winfileend
 character*6 winfileend
  
  !default fileending: non spin-polarized
  w2winfileend = ".w2win"
  winfileend = ".win"
  startmessage = "++ Preparing a non spin-polarized input file ++"
 
 !command line argument read-in
 iarg=iargc()
 argcount = 1
 if(iarg.ge.1) then
     do j=1,iarg
        call getarg(j,argdummy)
        if (argdummy(1:1).eq.'-') then
           if ((argdummy(2:3).eq.'up').or.(argdummy(2:3).eq.'dn')) then     
              !for spin-polarized calc. the fileendings have additional up/dn
              w2winfileend = ".w2win"//argdummy(2:3)
              winfileend = ".win"//argdummy(2:3)
              startmessage = "++ Preparing a spin-polarized input file:"//argdummy(2:3)//" ++"
           else
              write(*,*)"Error: unknown option"
              stop
           endif
        else
            if (argcount.eq.1) then
               read(argdummy,*)seedname
               argcount = argcount + 1
            else
               write(*,*)"Error: unknown option/input"
               stop
            endif
        endif
     enddo
  
  else
      write(*,*) 'Usage: write_win [-up/-dn] SEEDNAME'
      stop
  endif

   write(*,*)startmessage

  !read in number of Bloch bands
  open(unit=1,file=clearspace(seedname)//clearspace(w2winfileend),status='old')
  read(1,*)dummy
  read(1,"(2I3,A60)")nemin,nemax,dummy
  options%num_bands = nemax - nemin + 1
  read(1,"(2I3,A60)")tmpint,options%num_wann
  close(1)

 
 !read-in from seed.struct/seed.outputkgen
 open(unit=3,file=clearspace(seedname)//'.struct',status='old')
 open(unit=4,file=clearspace(seedname)//'.outputkgen',status='old')
 call countatoms(3,lattice)
 call readin_lattice(3,4,lattice)
 close(3)
 close(4)
 !write(*,*)lattice%transform_matrix

 !read-in from seed.klist
 open(unit=1,file=clearspace(seedname)//'.klist',status='old')
 !nkpoints = line_count(clearspace(seedname)//'.klist',1) - 1
 
 read(1,"(I10,3I10,I10,3F5.1,I10,A10,3I3)")tmp,vec,divisor,tmp2,tmp,dum2,mp_grid
 rewind(1)
 nkpoints = mp_grid(1)*mp_grid(2)*mp_grid(3)
 allocate(kpoints(nkpoints,3))

 do j1=1,nkpoints
   
   if (j1.eq.1) then
      read(1,"(I10,3I10,I10,3F5.1,I10,A10,3I3)")tmp,vec,divisor,tmp2,tmp,dum2,mp_grid
      !read(1,*)tmp,vec,divisor
   else
      read(1,*)tmp,vec
   endif
   do j2=1,3
      rvec(j2) = real(vec(j2))/real(divisor)
   enddo
   !rvec = matmul(lattice%transform_matrix,rvec)      
   !rvec = matmul(lattice%basis_real,rvec)    
   kpoints(j1,:) = rvec
  
 enddo
 write(*,*)"Found ",nkpoints," k-points from LAPW computation"
 close(1)

 !rvec(1)=0d0
 !rvec(2)=1d0
 !rvec(3)=0d0

  !write(*,*)lattice%basis_real
 !write(*,*)matmul(lattice%basis_real,rvec)
 !read-in k-path from seed.klist_band
 inquire(file=clearspace(seedname)//'.klist_band', exist=write_kpath) 
 if (write_kpath) then
    open(unit=4,file=clearspace(seedname)//'.klist_band',status='old')
    nkpoints_band = 0
    nlines=line_count(clearspace(seedname)//'.klist_band',4)
    
    call countkpoints_band(4,clearspace(seedname)//'.klist_band',kpoints_band,kpoints_band_spec)

    second_kpath_point = .true.
    counter = 0
    do j1=1,nlines-1
      read(4,"(A60)")dummy
      counter = counter + 1
   !   write(*,*)dummy
      if (clearspace(dummy(1:6)).ne."") then
         nkpoints_band = nkpoints_band + 1
         kpoints_band_spec(nkpoints_band) = dummy(1:6)
         write(*,*)dummy(11:60)
         read(dummy(11:60),"(4I5)")vec,divisor
         do j2=1,3
               rvec(j2) = real(vec(j2))/real(divisor)
         enddo
         !write(*,*)rvec
         !rvec = matmul(lattice%transform_matrix,rvec)
         !write(*,*)rvec
         !rvec = matmul(lattice%basis_real,rvec)
         !write(*,*)rvec
         kpoints_band(nkpoints_band,:) = rvec
         !write(*,*)kpoints_band(nkpoints_band,:)
         if ((counter.gt.1).and.(second_kpath_point.eqv..true.)) then
            nkpoints_firstpath = counter-1
            second_kpath_point = .false.
         endif
      endif
    enddo
    !write(*,*)kpath_density
    write(*,*)"Found ",nkpoints_band," k-points for bandstructure"
    close(4)
    
 endif
 !open(unit=4,file=clearspace(seedname)//'.klist_band',status='unknown')
 
 
 !write in seed.win
 open(unit=2,file=clearspace(seedname)//clearspace(winfileend),status='unknown')
 write(*,*)"Write in: ",clearspace(seedname)//clearspace(winfileend)
 
 !do j1=1,35
 !   read(2,*) dummy
 !   write(*,*) dummy
 !enddo
 write(2,*)"iprint = 3"
 write(2,"(A19,I3)")" num_bands       = ",options%num_bands
 write(2,"(A19,I3)")" num_wann        = ",options%num_wann
 write(2,*)"num_iter        = 1000"
 write(2,*)"num_print_cycles =100"
 write(2,*)""
 write(2,*)"!dis_froz_min     = 7."
 write(2,*)"!dis_froz_max     = 9."
 write(2,*)"dis_mix_ratio   = 0.5"
 write(2,*)"write_proj = .true."
 write(2,*)"write_xyz = .true."
 write(2,*)"translate_home_cell = .true."
 write(2,*)""
 write(2,*)"!SYSTEM"
 write(2,*)""
 write(2,*)"begin unit_cell_cart"
 write(2,*)lattice%unit
 !if (lattice%orthogonal) then
   write(2,"(3F11.7)")lattice%basis_real(1,:)
   write(2,"(3F11.7)")lattice%basis_real(2,:)
   write(2,"(3F11.7)")lattice%basis_real(3,:)
 !else
 !  write(*,*)"non-orthogonal lattices not implemented"
 !  stop
 !endif
 write(2,*)"end unit_cell_cart"

 write(2,*)"begin atoms_frac"
 atomcount=0
 do jelement=1,size(lattice%elements)
    do jmult=1,lattice%multiplicity(jelement)
       atomcount = atomcount + 1
       write(2,"(A2,3F8.5)")lattice%elements(jelement),lattice%atom_positions(atomcount,:)
    enddo
 enddo
 !Fe 0.00000 0.00000 0.00000 
 !Fe 0.50000 0.50000 0.50000
 !Sb 0.18800 0.35700 0.00000
 !Sb 0.81200 0.64300 0.00000
 !Sb 0.31200 0.85700 0.50000
 !Sb 0.68800 0.14300 0.50000
 write(2,*)"end atoms_frac"
 write(2,*)""
 write(2,*)"begin projections"
 !#V:l=0
 !#V:l=0
 !#V:l=0
 write(2,*)"end projections"
 write(2,*)""
 write(2,*)"begin kpoint_path"
 if (write_kpath) then
    do j1=1,nkpoints_band-1
      write(2,"(A3,3F6.2,A4,3F6.2)")clearspace(kpoints_band_spec(j1)),kpoints_band(j1,:),&
                                    clearspace(kpoints_band_spec(j1+1)),kpoints_band(j1+1,:)
    enddo
 endif
 !G 0.00  0.00  0.00    a 0.00  0.00  0.50
 !a 0.00  0.00  0.50    b 0.50  0.00  0.50
 !b 0.50  0.00  0.50    c 0.50  0.50  0.50
 !c 0.50  0.50  0.50    G 0.00  0.00  0.00
 !G 0.00  0.00  0.00    d 0.50  0.00  0.00
 !d 0.50  0.00  0.00    e 0.50  0.50  0.00
 !e 0.50  0.50  0.00    G 0.00  0.00  0.00
 write(2,*)"end kpoint_path"
 write(2,*) 
 write(2,*)"bands_plot = .false."
 write(2,"(A19,I3)")"bands_num_points  =",nkpoints_firstpath
 write(2,*)"wannier_plot = .false."
 write(2,*)"!restart = plot"
 write(2,*) 
 write(2,*)'! KPOINTS'  
 write(2,*)
 write(2,*)'mp_grid : ',mp_grid
 write(2,*)
 write(2,*)'begin kpoints'
 do j1=1,nkpoints
    write(2,"(3F13.9)") kpoints(j1,:)
 enddo
 write(2,*)'end kpoints'
 close(2)
 !close(12)
 deallocate(kpoints)
END PROGRAM



