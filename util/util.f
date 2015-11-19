MODULE util
!P.Wissgott 10/01/10
!*****************************************************    
   implicit none

   type structure
     logical orthogonal
     real*8 constants(3),angles(3)
     real*8 basis_real(3,3),transform_matrix(3,3)
     real*8, allocatable :: atom_positions(:,:)
     real*4, allocatable :: Z(:)
     integer, allocatable :: multiplicity(:)
     character*2, allocatable :: elements(:)
     character*1 specifier
     integer ::   atomcount,nop,elementcount
     character*4  :: unit
     real*8,allocatable :: sym(:,:,:)
   end type structure

   type optiontype
      integer :: num_bands
      integer :: num_wann
   end type optiontype
  
contains
!-----------------------------------------------------------------
  SUBROUTINE countatoms(unit,lattice)
  !counts the atom types and atoms in struct-file unit
  
  implicit none
  integer unit,mult,elementcount,atomcount,j
  logical readatoms
  character*60 dummy
  type(structure) lattice

  do j = 1,4
    read(unit,*)dummy
  enddo

  readatoms = .true.
  elementcount = 0
  atomcount = 0
  do while (readatoms)

    read(unit,"(A60)") dummy
    if (dummy(1:4).eq."ATOM") then
        elementcount = elementcount + 1
        atomcount = atomcount + 1
        read(unit,"(A15,I2)")dummy,mult
        do j = 1,mult-1
           atomcount = atomcount + 1
           read(unit,*)dummy
        enddo
    else   
       readatoms = .false.
    endif
    read(unit,"(A60)") dummy
    read(unit,"(A60)") dummy
    read(unit,"(A60)") dummy
    read(unit,"(A60)") dummy
   enddo
    
   lattice%elementcount = elementcount
   allocate(lattice%multiplicity(elementcount))
   allocate(lattice%elements(elementcount)) 
   allocate(lattice%Z(elementcount))
   allocate(lattice%atom_positions(atomcount,3))
   !write(*,*)elementcount,atomcount
   rewind(unit)

   END SUBROUTINE

!------------------------------------------------------------------

  SUBROUTINE readin_lattice(unit_struct,unit_outputkgen,lattice)

  implicit none
  
  integer jelement,jmult,jwrite,j,i !loop variables
  integer mult !multiplicity
  integer atomcount !counter for unit cell atoms
  integer nlines   !# lines in a file
  integer unit_struct,unit_outputkgen    !input units
  integer loccount
  character*60 dummy
  character*12 dum4
  type(structure) lattice
   
  !read(unit_struct,*)dummy
  read(unit_struct,*)dummy
  read(unit_struct,*)dummy
  read(dummy(1:1),*)lattice%specifier
  read(unit_struct,"(A23,A4)")dummy,lattice%unit
  if ((lattice%unit.ne."bohr").and.(lattice%unit.ne."ang"))  then
    write(*,*)"WARNING: did not recognize unit specifier in struct file: ",lattice%unit
    write(*,*)">>> set units to default value bohr"
    lattice%unit = "bohr"
  endif
  read(unit_struct,"(6F10.6)")lattice%constants,lattice%angles
  if ((lattice%angles(1).eq.90d0).and.(lattice%angles(2).eq.90d0).and.(lattice%angles(3).eq.90d0) ) then
     lattice%orthogonal = .true.
  else
     lattice%orthogonal = .false. 
  endif

    atomcount = 0 
  do jelement=1,size(lattice%elements)

     atomcount = atomcount + 1
   
     read(unit_struct,"(A12,F10.8,A3,F10.8,A3,F10.8)")dum4,lattice%atom_positions(atomcount,1),&
                                              dum4,lattice%atom_positions(atomcount,2),&
                                              dum4,lattice%atom_positions(atomcount,3)
     read(unit_struct,"(A15,I2)")dummy,mult
     lattice%multiplicity(jelement) = mult
  
     do jmult = 1,mult-1
            atomcount = atomcount + 1
            read(unit_struct,"(A60)")dummy
            read(dummy,"(A12,F10.8,A3,F10.8,A3,F10.8)")dum4,lattice%atom_positions(atomcount,1),&
                                                      dum4,lattice%atom_positions(atomcount,2),&
                                                      dum4,lattice%atom_positions(atomcount,3)
     enddo

     read(unit_struct,"(A60)") dummy
     lattice%elements(jelement) = dummy(1:2)
     read(dummy(56:59),*)lattice%Z(jelement)
     read(unit_struct,"(A60)") dummy
     read(unit_struct,"(A60)") dummy
     read(unit_struct,"(A60)") dummy

  enddo
  !read-in symmetry operations
  read(unit_struct,"(A60)")dummy
  read(dummy(1:7),"(I7)") lattice%nop
  write(*,*)"debug",lattice%nop
  allocate(lattice%sym(lattice%nop,3,3))
  do j=1,lattice%nop
      do i=1,3
         read(unit_struct,"(3F2.0,A12)")lattice%sym(j,i,:),dum4
      enddo
      read(unit_struct,*)
  enddo


  lattice%atomcount = atomcount
  
  !read in from seed.outputkgen
     read(unit_outputkgen,*)dummy
     read(unit_outputkgen,"(A7,3F10.6)")dummy,lattice%basis_real(1,:)
     read(unit_outputkgen,"(A7,3F10.6)")dummy,lattice%basis_real(2,:)
     read(unit_outputkgen,"(A7,3F10.6)")dummy,lattice%basis_real(3,:)
     lattice%transform_matrix(:,:) = 0d0
     do j=1,3
        lattice%transform_matrix(j,1) = lattice%basis_real(j,1)/lattice%constants(1)
        lattice%transform_matrix(j,2) = lattice%basis_real(j,2)/lattice%constants(2)
        lattice%transform_matrix(j,3) = lattice%basis_real(j,3)/lattice%constants(3)
     enddo
  write(*,*)"Found ",lattice%elementcount, " elements"
  write(*,*)"Found ",atomcount, " atoms:"
  loccount = 1
  do jwrite=1,size(lattice%elements)
     
     if (lattice%multiplicity(jwrite).gt.1) then
        write(*,"(I3,A1,I3,A1,A2,A7,I2,A3,F4.1)")loccount,'-',loccount+lattice%multiplicity(jwrite)-1,':',lattice%elements(jwrite),&
                                                 & ", mult=",lattice%multiplicity(jwrite),",Z=",lattice%Z(jwrite)
        loccount = loccount+lattice%multiplicity(jwrite)
     else
        write(*,"(I7,A1,A2,A7,I2,A3,F4.1)")loccount,':',lattice%elements(jwrite),", mult=",lattice%multiplicity(jwrite),&
                                        & ",Z=",lattice%Z(jwrite)
        loccount = loccount + 1
     endif

  enddo

  END SUBROUTINE readin_lattice
!-----------------------------------------------------------------

  SUBROUTINE count_kpoints(unit,filename,nkpoints)
  !counts the kpoints from .outputkgen file
  implicit none
  
  integer j,nkpoints,nlines
  integer unit
  integer markon, markoff
  character(len=*) filename
  character*70 :: dummy

  nlines=line_count(unit)
  
  do j=1,nlines
     read(unit,"(A70)")dummy
     !write(*,*)dummy(3:30)
     if (dummy(3:34).eq."internal and cartesian k-vectors") then
        !write(*,*)"readon"
         markon = j
     elseif (dummy(3:30).eq."NO. OF INEQUIVALENT K-POINTS") then
        !write(*,*)"readoff"
         markoff = j
     endif
  enddo
  nkpoints = markoff - markon - 1
  rewind(unit)

  END SUBROUTINE count_kpoints



!---------------------------------------------------------------





  SUBROUTINE countkpoints_band(unit,filename,kpoints,kpoint_spec)
  !counts the number of kpoints for bandstructure plot
  
  implicit none
  integer unit,nlines,nkpoints_band,j
  character*60 dummy
  character(len=*) filename
  real*8, allocatable :: kpoints(:,:)
  character*3, allocatable :: kpoint_spec(:)


  nlines=line_count(unit)
  nkpoints_band = 0
  do j=1,nlines
       read(4,"(A60)")dummy
       if (clearspace(dummy(1:6)).ne."") then
          nkpoints_band = nkpoints_band + 1
       endif
  enddo
  allocate(kpoints(nkpoints_band-1,3))
  allocate(kpoint_spec(nkpoints_band-1))
    
  
   
   rewind(unit)

   END SUBROUTINE
!-----------------------------------------------------------------
   FUNCTION clearspace(s1)  result (s2)
   !clears space after string and returns shorter string
   
   character(*) :: s1
   character(char_count(s1)) :: s2
   integer :: i, n

   n = 0
   do i = 1,len_trim(s1)
      if (s1(i:i) == ' ') cycle
         n = n+1
         s2(n:n) = s1(i:i)
   end do

   END FUNCTION
!****************************************************

   PURE FUNCTION char_count(s)  result (n)
   !counts characters of string

   character(*),intent(in) :: s
   integer :: n, i
   
   n = 0
   do i = 1,len_trim(s)
      if (s(i:i) == ' ') cycle
        n = n+1
   end do
   
   END FUNCTION
!--------------------------------------------------
   FUNCTION get_string(n)  result (s)
   !convert integer to string

   integer,intent(in) :: n
   character*5  :: s
   
   
   write(s,"(I5)")n
   
   END FUNCTION

   integer function line_count(fid)  
   !>Returns the number of lines in a file.
   !>\param fname Name of the file
   !>\return Number of lines in the file.
	
   implicit none

   !input parameters
   !character(len=*) fname   !filename of the file to count
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


END MODULE util
