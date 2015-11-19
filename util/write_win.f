PROGRAM write_win
 !program prepares input case.win for wannier90 
 ! P.Wissgott 10/01/10
 
 use util 

 implicit none
 character*50 seedname
 integer nkpoints,j1,j2,tmp
 integer divisor,mp_grid(3),mp_grid_fine(3)!k-mesh specifiers
 integer iarg,argcount !command line input argument counter
 integer vec(3)
 integer jelement,jmult,jwrite !loop variables
 integer atomcount !counter for unit cell atoms
 integer nlines   !# lines in a file
 integer nkpoints_band !#kpoints for bandstructure plots
 integer nemin,nemax,tmpint !read in dummys
 integer counter,j
 integer nkpoints_firstpath  !number of kpoints along the first part of the kpath
 integer unitkgen,unitstruct, unitklist, unitw2win,unitkfine
 integer kcounter
 real*8, allocatable :: kpoints_band(:,:),kpoints(:,:),kpoints_fine(:,:) !k-mesh points 
 real*8 :: rvec(3)
 character*3, allocatable :: kpoints_band_spec(:),kpoints_fine_spec(:)
 real*8 :: tmp2(3), a1,a2,a3
 character*60 dummy
 character*10 dum2
 character*4 dum3!,unit
 character*70 argdummy,startmessage

 logical write_kpath,second_kpath_point
 logical :: kread,lattconv

 logical :: update,kfine
 real*8  :: dis_froz_min,dis_froz_max
 integer :: nkpoints_fine
 character*10 :: restart

 type(structure) lattice
 type(optiontype) options
 
 character*8 w2winfileend
 character*6 winfileend
  
  !setunit for file io
  unitstruct = 1
  unitkgen = 3
  unitw2win = 4
  unitklist = 100
  unitkfine = 11

  update = .false.
  lattconv = .false.

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
           !elseif ((argdummy(2:5).eq.'soup').or.(argdummy(2:5).eq.'sodn')) then     
           !   !for spin-polarized calc. the fileendings have additional up/dn
           !   w2winfileend = ".w2win"//argdummy(4:5)
           !   winfileend = ".win"//argdummy(2:3)
           !   startmessage = "++ Preparing a spin-polarized input file for SO calc.:"//argdummy(4:5)//" ++"
           elseif (argdummy(2:7).eq.'renew') then
               update = .true.

           elseif (argdummy(2:2).eq.'l') then
               lattconv = .true.
              
           elseif (argdummy(2:2).eq.'h') then
               write(*,*) 'prepares input file case.win[up/dn] for wannier90'
               write(*,*) 'Usage: write_w2win [-up/-dn] case'
               stop 
           else
              write(*,*)"Usage: write_win [-up/-dn] SEEDNAME"
              stop
           endif
        else
            if (argcount.eq.1) then
               read(argdummy,*)seedname
               argcount = argcount + 1
            else
               write(*,*)"Usage: write_win [-up/-dn] SEEDNAME"
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
  open(unit=unitw2win,file=clearspace(seedname)//clearspace(w2winfileend),status='old')
  read(unitw2win,*)dummy
  read(unitw2win,"(2I4,A60)")nemin,nemax,dummy
  options%num_bands = nemax - nemin + 1
  read(unitw2win,"(2I3,A60)")tmpint,options%num_wann
  close(unitw2win)

 
 !read-in from seed.struct/seed.outputkgen
 open(unit=unitstruct,file=clearspace(seedname)//'.struct',status='old')
 open(unit=unitkgen,file=clearspace(seedname)//'.outputkgen',status='old')
 call countatoms(unitstruct,lattice)
 call readin_lattice(unitstruct,unitkgen,lattice)
!  call count_kpoints(unitkgen,clearspace(seedname)//'.outputkgen',nkpoints)
 
 
 !readin kpoints
 nlines=line_count(unitkgen)
 kread = .false.
 kcounter = 1

!  do j=1,nlines
!      read(unitkgen,"(A60)")dummy
!      if (dummy(3:34).eq."internal and cartesian k-vectors") then
!          kread = .true.
!          cycle
!      elseif (dummy(3:30).eq."NO. OF INEQUIVALENT K-POINTS") then
!          kread = .false.
!          cycle
!      elseif (.not.kread) then
!          cycle
!      endif
!      if (kread) then
!         read(dummy(1:36),"(3F12.10)")kpoints(kcounter,1:3) 
!         kcounter = kcounter + 1
!      endif
!  enddo
 
 close(unitstruct)
 close(4)
 
 !read-in from seed.klist
 open(unit=unitklist,file=clearspace(seedname)//'.klist_w90',status='old')
 nkpoints = line_count(unitklist) -2 !take care
 write(*,*)"Found ",nkpoints," k-points from LAPW computation"
 allocate(kpoints(nkpoints,3))
 read(unitklist,"(I10,3I10,I10,3F5.1,I10,A10,3I3)")tmp,vec,divisor,tmp2,tmp,dum2,mp_grid
 do j2=1,3
     rvec(j2) = real(vec(j2))/real(divisor)
 enddo
 if ((lattice%specifier.ne."R").and.(lattice%specifier.ne."H").and.(lattice%specifier.ne."C").and.(.not.lattconv)) then
   rvec = matmul(lattice%transform_matrix,rvec)
 endif
 kpoints(1,:) = rvec
 do j1=2,nkpoints
    read(unitklist,"(A60)")dummy
      !write(*,*) dummy(11:52)
    read(dummy(11:52),*)vec,divisor
      !write(*,*)vec,divisor
    do j2=1,3
        rvec(j2) = real(vec(j2))/real(divisor)
    enddo
    if ((lattice%specifier.ne."R").and.(lattice%specifier.ne."H").and.(lattice%specifier.ne."C").and.(.not.lattconv)) then
         rvec = matmul(lattice%transform_matrix,rvec)
    endif
    !rvec = matmul(lattice%transform_matrix,rvec)
    kpoints(j1,:) = rvec
  enddo
  close(unitklist)
 

 inquire(file=clearspace(seedname)//'.klist',exist=kfine)
 kfine = .false. !turn off additional win part
 if (kfine) then
    open(unit=unitkfine,file=clearspace(seedname)//'.klist',status='old') 
    nlines = line_count(unitkfine)
    !call countkpoints_band(unitkfine,clearspace(seedname)//'.klist_fine',kpoints_fine,kpoints_fine_spec)
    read(unitkfine,"(I10,3I10,I10,3F5.1,I10,A10,3I3)")tmp,vec,divisor,tmp2,tmp,dum2,mp_grid_fine
    write(*,*)tmp,vec,divisor,tmp2,tmp,dum2,mp_grid_fine
    nkpoints_fine = nlines-2!take care!mp_grid_fine(1)*mp_grid_fine(2)*mp_grid_fine(3)
    mp_grid_fine = 1
    mp_grid_fine(1) = nkpoints_fine
    allocate(kpoints_fine(nkpoints_fine,3))
    write(*,*)"Found ",nkpoints_fine," fine k-points"
!     write(*,*)mp_grid_fine
!     write(*,*)"rvec1",vec,divisor
    do j2=1,3
        rvec(j2) = real(vec(j2))/real(divisor)
    enddo
!     write(*,*)"rvec2",rvec
    !stop
    if ((lattice%specifier.ne."R").and.(lattice%specifier.ne."H").and.(lattice%specifier.ne."C").and.(.not.lattconv)) then
       rvec = matmul(lattice%transform_matrix,rvec)
    endif
!     write(*,*)"rvec3",rvec
    kpoints_fine(1,:) = rvec
    do j1=2,nkpoints_fine
      read(unitkfine,"(A60)")dummy
      !write(*,*) dummy(11:52)
      read(dummy(11:52),*)vec,divisor
!       write(*,*)vec,divisor
      do j2=1,3
          rvec(j2) = real(vec(j2))/real(divisor)
      enddo
      if ((lattice%specifier.ne."R").and.(lattice%specifier.ne."H").and.(lattice%specifier.ne."C").and.(.not.lattconv)) then
            rvec = matmul(lattice%transform_matrix,rvec)
      endif
      kpoints_fine(j1,:) = rvec
    enddo
    close(unitkfine)
!     do j1=1,nkpoints_fine
!         write(*,*)"write",kpoints_fine(j1,:)
!     enddo
 endif
 
 
 
 
 !read-in k-path from seed.klist_band
 inquire(file=clearspace(seedname)//'.klist_band', exist=write_kpath) 
 if (write_kpath) then
    open(unit=4,file=clearspace(seedname)//'.klist_band',status='old')
    nkpoints_band = 0
    nlines=line_count(4)
    call countkpoints_band(4,clearspace(seedname)//'.klist_band',kpoints_band,kpoints_band_spec)
    second_kpath_point = .true.
    counter = 0
    do j1=1,nlines-1
      read(4,"(A60)")dummy
      counter = counter + 1
      if (clearspace(dummy(1:6)).ne."") then
         nkpoints_band = nkpoints_band + 1
         kpoints_band_spec(nkpoints_band) = dummy(1:6)
         read(dummy(11:60),"(4I5)")vec,divisor
!          write(*,*)"see",dummy
         do j2=1,3
               rvec(j2) = real(vec(j2))/real(divisor)
         enddo
!          write(*,*)"see1",rvec
         if ((lattice%specifier.ne."R").and.(lattice%specifier.ne."H").and.(lattice%specifier.ne."C").and.(.not.lattconv)) then
            rvec = matmul(lattice%transform_matrix,rvec)
         endif
         kpoints_band(nkpoints_band,:) = rvec
!          write(*,*)"see2",rvec
         if ((counter.gt.1).and.(second_kpath_point.eqv..true.)) then
            nkpoints_firstpath = counter-1
            second_kpath_point = .false.
         endif
      endif
    enddo
    write(*,*)"Found ",nkpoints_band," k-points for bandstructure"
    close(4)
    
 endif
 

 if (update) then
    dis_froz_min = 0d0
    dis_froz_max = 0d0
    restart = ""
    dummy = ""
    open(unit=2,file=clearspace(seedname)//clearspace(winfileend),status='old')
    nlines=line_count(2)
    write(*,*)nlines
    do j1=1,nlines
      read(2,"(A60)") dummy
      !write(*,*)dummy
      if (dummy(2:13).eq."dis_froz_min") then
         dummy = trim(dummy)
         read(dummy,*)dum2,dum3,dis_froz_min
      endif
      if (dummy(2:13).eq."dis_froz_max") then
         dummy = trim(dummy)
         read(dummy,*)dum2,dum3,dis_froz_max
      endif
      if (dummy(2:9).eq."restart") then
         dummy = trim(dummy)
         read(dummy,*)dum2,dum3,restart
      endif
    enddo
    close(2)
 endif
 !write in seed.win
 open(unit=2,file=clearspace(seedname)//clearspace(winfileend),status='unknown')
 write(*,*)"Write in: ",clearspace(seedname)//clearspace(winfileend)
 
 write(2,*)"iprint = 3"
 write(2,"(A19,I3)")" num_bands       = ",options%num_bands
 write(2,"(A19,I3)")" num_wann        = ",options%num_wann
 write(2,*)"num_iter        = 1000"
 write(2,*)"num_print_cycles =100"
 write(2,*)"!conv_window = 20"
 write(2,*)"!conv_tol = 0.0001"
 write(2,*)"!conv_noise_amp = 1"
 write(2,*)"!conv_noise_num = 3"
 write(2,*)""
 if (update) then
   write(*,"(A20,F7.3)")" dis_froz_min     = ",dis_froz_min
   write(*,"(A20,F7.3)")" dis_froz_max     = ",dis_froz_max
   write(2,"(A20,F7.3)")" dis_froz_min     = ",dis_froz_min
   write(2,"(A20,F7.3)")" dis_froz_max     = ",dis_froz_max
 else
   write(2,*)"!dis_froz_min     = 7."
   write(2,*)"!dis_froz_max     = 9."
 endif
 write(2,*)"dis_mix_ratio   = 0.5"
 write(2,*)"write_proj = .true."
 write(2,*)"write_xyz = .true."
 write(2,*)"translate_home_cell = .true."
 write(2,*)""
 write(2,*)"!SYSTEM"
 write(2,*)""
 write(2,*)"begin unit_cell_cart"
 write(2,*)lattice%unit
 write(2,"(3F12.7)")lattice%basis_real(1,:)
 write(2,"(3F12.7)")lattice%basis_real(2,:)
 write(2,"(3F12.7)")lattice%basis_real(3,:)
 write(2,*)"end unit_cell_cart"
 write(2,*)""
 write(2,*)"begin atoms_cart"
 atomcount=0
 do jelement=1,size(lattice%elements)
    do jmult=1,lattice%multiplicity(jelement)
       atomcount = atomcount + 1
       do j1=1,3
          rvec(j1) = lattice%constants(j1)*lattice%atom_positions(atomcount,j1)
       enddo
       write(2,"(A2,3F9.5)")lattice%elements(jelement),rvec(:)*0.529177249d0!lattice%atom_positions(atomcount,:)
    enddo
 enddo
 write(2,*)"end atoms_cart"
 write(2,*)""
 write(2,*)"begin projections"
 write(2,*)"end projections"
 write(2,*)""
 write(2,*)"begin kpoint_path"
 if (write_kpath) then
    do j1=1,nkpoints_band-1
      write(2,"(A3,3F6.2,A4,3F6.2)")clearspace(kpoints_band_spec(j1)),kpoints_band(j1,:),&
                                    clearspace(kpoints_band_spec(j1+1)),kpoints_band(j1+1,:)
    enddo
 endif
 write(2,*)"end kpoint_path"
 write(2,*) 
 if (write_kpath) then
     write(2,*)"bands_plot = .true."
     write(2,"(A20,I3)")" bands_num_points  =",nkpoints_firstpath
     write(2,*)"!bands_plot_mode = cut"
     write(2,*)"!bands_plot_project = 1"
 else
     write(2,*)"bands_plot = .false."
 endif
 write(2,*)"wannier_plot = .false."
 if ((update).and.(restart.ne."")) then
    write(2,"(A11,A10)")" restart = ",trim(restart)
 else
    write(2,*)"!restart = plot"
 endif
 write(2,*)"hr_plot = .true."
 write(2,*)"!dist_cutoff = 10"
 write(2,*) 
 write(2,*)"kmesh_tol = 0.0001"
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
 if (kfine) then
    write(2,*)
    write(2,*)'mp_fine : ',mp_grid_fine
    write(2,*)
    write(2,*)'begin kpfine'
    do j1=1,nkpoints_fine
       write(2,"(3F13.9)") kpoints_fine(j1,:)
    enddo
    write(2,*)'end kpfine'
 endif
 close(2)
 deallocate(kpoints)
END PROGRAM



