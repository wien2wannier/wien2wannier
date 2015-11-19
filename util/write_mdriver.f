PROGRAM write_mdriver
!Program writes matlab-driver file for visualization
! P.Wissgott 10/01/10

use util

implicit none

character*50   :: seedname
character*100  :: seedpath
character*60   :: dummy
integer        :: iarg
integer        :: unitwin,unitd,unitstruct,unitin7,unitoutputkgen
integer        :: num_wann
integer        :: jgrid,jwrite
integer        :: divisor
real*8         :: Zdivisor
real*8         :: basis_cart(3,3)
real*8         :: grid_origin(3),grid_vec(3,3)
integer        :: tmp(3)
character*20   :: frmt

 
type(structure) lattice

  !command line argument read-in
  iarg=iargc()
  if(iarg.eq.1) then
     call getarg(1,seedname)
     seedpath = ''
  elseif (iarg.eq.2) then
     call getarg(1,seedname)
     call getarg(2,seedpath)
  else
      write(*,*) 'Error: at least commandline argument "case" has to be given.'
      stop
  endif

  unitstruct = 11
  unitwin = 12
  unitin7 = 13
  unitd = 14
  unitoutputkgen = 15

  !read-in from struct file
  !remark: we take the conventional basis from struct file!
  open(unit=unitstruct,file=trim(seedpath)//trim(seedname)//".struct",status="old")
  open(unit=unitoutputkgen,file=trim(seedpath)//trim(seedname)//".outputkgen",status="old")
  call countatoms(unitstruct,lattice)
  call readin_lattice(unitstruct,unitoutputkgen,lattice)
  close(unitstruct)  
  basis_cart(:,:) = 0d0
  basis_cart(1,1) = lattice%constants(1)
  basis_cart(2,2) = lattice%constants(2)
  basis_cart(3,3) = lattice%constants(3)
  !write(*,*)"debug",lattice%constants,lattice%unit
 
 
  !convert to Angstrom
  if (lattice%unit.eq."bohr") then
      basis_cart = 0.5292*basis_cart
  endif
  write(*,*)"Use cart. basis from ",trim(seedname)//".win (in Anstrom)"
  write(*,"(3F8.4)")basis_cart(:,1)
  write(*,"(3F8.4)")basis_cart(:,2)
  write(*,"(3F8.4)")basis_cart(:,3)
  !close(unit=unitwin)

  !read-in from in7 file(lapw7)
  open(unit=unitin7,file=trim(seedpath)//trim(seedname)//".wplotin",status="old")
  read(unitin7,*)dummy
  read(unitin7,"(I2,I3,I3,I2)")tmp,divisor
  grid_origin = real(tmp)/real(divisor)
  write(*,*)"Use spatial box from ",trim(seedname)//".wplotin (in cart. basis vectors)"
  write(*,"(A14,3F7.2)")"origin: ",grid_origin
  
  do jgrid=1,3
     read(unitin7,"(I2,I3,I3,I2)")tmp,divisor
     grid_vec(jgrid,:) = real(tmp)/real(divisor)
     write(*,"(I1,A13,3F7.2)")jgrid,"-D-boundary: ",grid_vec(jgrid,:) 
  enddo
  close(unitin7)
  
  !write driver file
  open(unit=unitd,file=trim(seedname)//"_mdriver.m",status="unknown") 
  write(unitd,*)"function out="//trim(seedname)//"_mdriver(flag,varargin)"
  write(unitd,*)"                                                               "
  write(unitd,"(A56)")"  switch (flag)                                          "
  write(unitd,"(A56)")"     case 'material'                                     "
  frmt = "(A12,A"//clearspace(get_string(char_count(seedname)))//",A2)"
  !write(*,*) frmt
  write(unitd,trim(frmt))"      out='",trim(seedname),"';"
  write(unitd,"(A56)")"                                                         "
  write(unitd,"(A56)")"     case 'path'                                         "
  if (char_count(seedpath).gt.0) then
    frmt = "(A14,A"//clearspace(get_string(char_count(seedpath)))//",A2)"
    write(unitd,trim(frmt))"         out='",trim(seedpath),"';"
  else
    write(unitd,"(A14)")"       out='';"
  endif
  write(unitd,"(A56)")"                                                         "
  write(unitd,"(A56)")"     case 'num_wann'                                     "
  write(unitd,"(A11,I2,A1)")"      out=",num_wann,";"
  write(unitd,"(A56)")"                                                         "
  write(unitd,"(A56)")"     case 'basis'                                        "
  write(unitd,"(A12,3F11.7,A1)")"      out=[",basis_cart(1,:),";"
  write(unitd,"(A12,3F11.7,A1)")"           ",basis_cart(2,:),";"
  write(unitd,"(A12,3F11.7,A2)")"           ",basis_cart(3,:),"];"
  write(unitd,"(A56)")"                                                         "
  write(unitd,"(A56)")"     case 'atom_positions'                               "
  write(unitd,"(A12,3F11.7,A1)")"      out=[",lattice%atom_positions(1,:),";"
  do jwrite=2,lattice%atomcount-1
     write(unitd,"(A12,3F11.7,A1)")"          ",lattice%atom_positions(jwrite,:),";"
  enddo
  write(unitd,"(A12,3F11.7,A2)")"           ",lattice%atom_positions(lattice%atomcount,:),"];"
  write(unitd,"(A56)")"                                                         "
  write(unitd,"(A56)")"     case 'atom_radii'                                   "
  Zdivisor = 70d0
  write(unitd,"(A12,F11.7,A1)")"      out=[",lattice%Z(1)/Zdivisor,";"
  do jwrite=2,size(lattice%elements)-1
     write(unitd,"(A12,F11.7,A1)")"          ",lattice%Z(jwrite)/Zdivisor,";"
  enddo
  write(unitd,"(A12,F11.7,A2)")"           ",lattice%Z(size(lattice%elements))/Zdivisor,"];"
  write(unitd,"(A56)")"                                                         "
  write(unitd,"(A56)")"     case 'multiplicities'                               "
  write(unitd,"(A12,I2,A1)")"      out=[",lattice%multiplicity(1),";"
  do jwrite=2,size(lattice%elements)-1
     write(unitd,"(A12,I2,A1)")"          ",lattice%multiplicity(jwrite),";"
  enddo
  write(unitd,"(A12,I2,A2)")"           ",lattice%multiplicity(size(lattice%elements)),"];"
  write(unitd,"(A56)")"                                                         "
  write(unitd,"(A56)")"     case 'grid_origin'                                  "
  write(unitd,"(A12,3F11.7,A2)")"      out=[",grid_origin,"];"
  write(unitd,"(A56)")"                                                         "
  write(unitd,"(A56)")"     case 'grid_vectors'                                 "
  write(unitd,"(A12,3F11.7,A1)")"      out=[",grid_vec(1,:),";"
  write(unitd,"(A12,3F11.7,A1)")"           ",grid_vec(2,:),";"
  write(unitd,"(A12,3F11.7,A2)")"           ",grid_vec(3,:),"];"
  write(unitd,"(A56)")"                                                         "
  write(unitd,"(A56)")"     %case 'axis'%in Cartesian coordinates of grid       "
  write(unitd,"(A56)")"     %out=[0 1 0 1 0 1];%manual axis definition          "
  write(unitd,"(A56)")"                                                         "
  write(unitd,"(A56)")"     case 'plotmode'                                     "
  write(unitd,"(A19)")"       out='plain';"
  write(unitd,"(A56)")"                                                         "
  write(unitd,"(A56)")"     case 'printmode'                                     "
  write(unitd,"(A19)")"       out='off';"
  write(unitd,"(A56)")"  end                                                    "
  close(unitd)
  write(*,*)"Wrote to "//trim(seedname)//"_mdriver.m"

END PROGRAM write_mdriver