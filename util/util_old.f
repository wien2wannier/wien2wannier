MODULE util
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
     integer ::   atomcount
     character*4  :: unit
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
        read(unit,"(A16,I2)")dummy,mult
        !write(*,*)dummy,mult
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
   !write(*,"(A60)")dummy
   enddo
    
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
  
  integer jelement,jmult,jwrite !loop variables
  integer mult !multiplicity
  integer atomcount !counter for unit cell atoms
  integer nlines   !# lines in a file
  integer unit_struct,unit_outputkgen    !input units
  character*60 dummy
  character*12 dum4
  type(structure) lattice
   
  !read(unit_struct,*)dummy
  read(unit_struct,*)dummy
  read(unit_struct,*)dummy
  read(unit_struct,"(A23,A4)")dummy,lattice%unit
  read(unit_struct,"(6F10.6)")lattice%constants,lattice%angles
  if ((lattice%angles(1).eq.90d0).and.(lattice%angles(2).eq.90d0).and.(lattice%angles(3).eq.90d0) ) then
     lattice%orthogonal = .true.
  else
     lattice%orthogonal = .false. 
  endif

    atomcount = 0 
  do jelement=1,size(lattice%elements)

     !read(unit_struct,"(A60)") dummy
    
     atomcount = atomcount + 1
   
     read(unit_struct,"(A12,F10.8,A3,F10.8,A3,F10.8)")dum4,lattice%atom_positions(atomcount,1),&
                                              dum4,lattice%atom_positions(atomcount,2),&
                                              dum4,lattice%atom_positions(atomcount,3)
     read(unit_struct,"(A16,I2)")dummy,mult
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
  lattice%atomcount = atomcount
  
  !read in from seed.outputkgen
  !if (lattice%orthogonal) then
     read(unit_outputkgen,*)dummy
     read(unit_outputkgen,"(A8,3F10.6)")dummy,lattice%basis_real(1,:)
     read(unit_outputkgen,"(A8,3F10.6)")dummy,lattice%basis_real(2,:)
     read(unit_outputkgen,"(A8,3F10.6)")dummy,lattice%basis_real(3,:)
     lattice%transform_matrix(:,:) = 0d0
     lattice%transform_matrix(1,1) = 1d0/lattice%constants(1)
     lattice%transform_matrix(2,2) = 1d0/lattice%constants(2)
     lattice%transform_matrix(3,3) = 1d0/lattice%constants(3)
  !else
  ! write(*,*)"non-orthogonal lattices not implemented"
  ! stop
 !endif

  write(*,*)"Found ",atomcount, " atoms:"
  do jwrite=1,size(lattice%elements)

     write(*,"(2A2,I2,A3,F4.1)")lattice%elements(jwrite),": ",lattice%multiplicity(jwrite),",Z=",lattice%Z(jwrite)

  enddo

  END SUBROUTINE readin_lattice
!-----------------------------------------------------------------
  SUBROUTINE countkpoints_band(unit,filename,kpoints,kpoint_spec)
  !counts the number of kpoints for bandstructure plot
  
  implicit none
  integer unit,nlines,nkpoints_band,j
  character*60 dummy
  character(len=*) filename
  real*8, allocatable :: kpoints(:,:)
  character*3, allocatable :: kpoint_spec(:)


  nlines=line_count(filename,unit)
  nkpoints_band = 0
  do j=1,nlines
       read(4,"(A60)")dummy
       if (clearspace(dummy(1:6)).ne."") then
          nkpoints_band = nkpoints_band + 1
          !write(*,*)dummy(1:6) 
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
   !write(*,*)n
   !write(*,*)s
   
   END FUNCTION

   integer function line_count(fname,fid)  
   !>Returns the number of lines in a file.
   !>\param fname Name of the file
   !>\return Number of lines in the file.
	
   implicit none

   !input parameters
   character(len=*) fname   !filename of the file to count
   integer fid 

   !local parameters
   character*20 dummy
   integer ioStat
   logical ioEndFlag

   ioEndFlag = .false.
   line_count = 0
		
   do while (.not. ioEndFlag )
      read(fid,"(A20)", iostat=ioStat ) dummy
!write(*,*)dummy
      if( iostat .eq. 0) line_count = line_count + 1
!write(*,*)line_count
      if( iostat < 0 ) ioEndFlag = .TRUE. 
   end do
   rewind(fid)

   END FUNCTION line_count
END MODULE util