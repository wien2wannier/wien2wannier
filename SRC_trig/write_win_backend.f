!!! wien2wannier/SRC_trig/write_win.f
!!!
!!!    Prepares input case.win for Wannier90
!!!
!!! Copyright 2009-2012 Philipp Wissgott
!!!           2013-2014 Elias Assmann
!!!
!!! $Id: write_win_backend.f 289 2014-10-10 12:18:40Z assmann $

program write_win
  use structmod, only: struct, struct_read
  use const,     only: BUFSZ, DPk
  use util,      only: lowercase, newunit
  use clio,      only: fetcharg, argstr
  use kpoints,   only: get_kmesh_band, get_kmesh_klist

  implicit none

!!!----------- Configurable parameters -----------
  !! Initial size and size increment of keys, vals.
  integer, parameter :: NKEY_FIRST=50, NKEY_INC=50

!!!----------- Variables               -----------
  integer :: unit_inwf=13

  integer :: nemin, nemax, num_wann, num_bands, mp_grid(3)
  integer :: c, i, j, k, nlin, nfill, tmpint, iarg

  real(DPk), allocatable :: kpath(:,:), kmesh(:,:)
  character, allocatable :: knames(:)
  integer,   allocatable :: centers(:)

  type(argstr) :: inwffile, structfile, klistfile, bandfile

  character(len=BUFSZ) :: line, key, filler, comment

  character(len=BUFSZ), allocatable :: keys(:), vals(:)
  logical,              allocatable :: keys_done(:)
  integer                           :: nkeys

  logical :: atoms_done = .false., uc_done = .false., proj_done = .false.
  logical :: kmesh_done = .false., kpath_done = .false.
  logical :: write_kpath, guiding_centres=.true.

  type(struct) :: stru

!!!----------- Code                    -----------
  call fetcharg(1, inwffile,   "failed to get `inwf' argument")
  call fetcharg(2, structfile, "failed to get `struct' argument")
  call fetcharg(3, klistfile,  "failed to get `klist' argument")
  call fetcharg(4, bandfile,   "failed to get `klist_band' argument")

!!! Read ‘inwf’ file for num_bands, num_wann
  open(newunit(unit_inwf), FILE=inwffile%s, STATUS='old')
  read(unit_inwf, *)
  read(unit_inwf, *) nemin, nemax
  num_bands = nemax - nemin + 1

  read(unit_inwf, *) tmpint, num_wann
  allocate(centers(num_wann))
  proj: do i = 1,num_wann
     read(unit_inwf, *) nlin
     ylm: do j = 1,nlin
        read(unit_inwf, *) c
        if (j==1) then
           centers(i) = c
        elseif (centers(i) /= c) then
           ! disable “guiding centres” if a projection is not uniquely
           ! centered on one atom
           guiding_centres = .false.
           exit proj
        end if
     end do ylm
  end do proj
  close(unit_inwf)

!!! Read ‘struct’
  call struct_read(structfile%s, stru)

!!! Read ‘klist’ for MP k-mesh
  call get_kmesh_klist(klistfile%s, kmesh, stru, mp_grid)

!!! Read parameters from command line
  call allocate_keyval(NKEY_FIRST)

  !! num_wann, num_bands, and mp_grid must match with ‘inwf’ and
  !! ‘klist’
        keys(1)        = 'num_wann'
  write(vals(1), '(I0)')  num_wann
        keys(2)        = 'num_bands'
  write(vals(2), '(I0)')  num_bands

  nkeys = 2
  optarg: do iarg = 5, command_argument_count(), 2
     nkeys = nkeys+1
     if (nkeys > size(keys)) call realloc_keyval(NKEY_INC)

     call fetcharg(iarg,   keys(nkeys))
     call fetcharg(iarg+1, vals(nkeys))
  end do optarg

  i = findkey('guiding_centres')
  if (i /= 0) read(vals(i), *) guiding_centres

!!! Read ‘klist_band’ for BZ path
  inquire(FILE=bandfile%s, EXIST=write_kpath)

  if (write_kpath) then
     call get_kmesh_band(bandfile%s, kpath, stru, knames)
  end if

!!! Reproduce template, putting in the values we got
  allocate(keys_done(nkeys))
  keys_done(:) = .false.

  template: do
     read (*, '(A)', END=101) line

     block: if (index(adjustl(lowercase(line)), 'begin') == 1) then
        i = index(line, 'begin')

        if      (index(adjustl(lowercase(line(i+6:))), 'atoms_cart') ==1) then
           call print_atoms(cart=.true.);  atoms_done = .true.
        else if (index(adjustl(lowercase(line(i+6:))), 'atoms_frac') ==1) then
           call print_atoms(cart=.false.); atoms_done = .true.
        else if (index(adjustl(lowercase(line(i+6:))), 'unit_cell')  ==1) then
           call print_uc();    uc_done    = .true.
        else if (index(adjustl(lowercase(line(i+6:))), 'projections')==1) then
           call print_proj();  proj_done  = .true.
        else if (index(adjustl(lowercase(line(i+6:))), 'kpoints')    ==1) then
           call print_kmesh(); kmesh_done = .true.
        else if (index(adjustl(lowercase(line(i+6:))), 'kpoint_path')==1) then
           call print_kpath(); kpath_done = .true.
        else
           call print_block()
        end if
        cycle template
     end if block

     j = scan(line, '#!')
     if (j==0) then
        j = len_trim(line)+1
     elseif (j>1) then
        j = verify(line(1:j-1), ' 	', BACK=.true.)+1
     end if
     comment = line(j:)

     i = verify(line(1:j-1), ' 	')
     empty: if (i==0) then
        print '(A)', trim(line)
        cycle template
     end if empty

     k = scan(line(i:j-1), ' 	=:')
     key    = line(i:i+k-2)
     j = verify(line(i+k-1:), ' 	=:')
     filler = line(i+k-1 : i+k+j-3)
     nfill  = j-1

     k = findkey(key)
     if (k==0) then
        print '(A)', trim(line)
     else
        j = nfill-len_trim(vals(k))
        print'(5A)', line(1:i-1), trim(key), filler(1:nfill), &
             &       trim(vals(k)), trim(comment)
        keys_done(k) = .true.
     end if
  end do template
101 continue

!!! Add any parameters that were not in the template
  if (.not. all(keys_done)) print*

  rest: do i=1,nkeys
     if (.not. keys_done(i)) &
          call printkey(i)
  end do rest

  if (.not.    uc_done) call print_uc   (append=.true.)
  if (.not. atoms_done) call print_atoms(append=.true.)
  if (.not. kpath_done) call print_kpath(append=.true.)
  if (.not.  proj_done) call print_proj (append=.true.)
  if (.not. kmesh_done) call print_kmesh(append=.true.)

!!!----------- Done.                   -----------

contains

!!!----------- Procedures for blocks we care about
  subroutine print_uc(append)
    logical, intent(in), optional :: append
    logical                       :: a
    character(len=*), parameter   :: head = 'begin unit_cell_cart'
    character(len=*), parameter   :: tail = 'end unit_cell_cart'

    a = .false.
    if (present(append)) a = append

    appending: if (a) then
       print*
       print '(A)', head
       print '(A)', '  Bohr'
    else
       print '(A)', trim(line)
       
       read(*, '(A)') line
       unit: if (index(adjustl(lowercase(line)), 'bohr') == 1) then
          print '(A)', trim(line)
       else
          print '(A)', '  Bohr'
       end if unit
    end if appending

    print '(3F12.7)', transpose(stru%brlat)

    if (a) then
       print '(A)', tail
    else
       call skip_to_end()
    end if
  end subroutine print_uc

  subroutine print_proj(append)
    logical, intent(in), optional :: append
    logical                       :: a
    character(len=*), parameter   :: head = 'begin projections'
    character(len=*), parameter   :: tail = 'end projections'
    integer                       :: i

    a = .false.
    if (present(append)) a = append

    if (a) then
       print*
       print '("!!!! Dummy `projections'' block for guiding centres !!!")'

       if (guiding_centres) then
          i = findkey('guiding_centres') ! must not print it twice
          if (i==0) print '(A)', 'guiding_centres = .true.'
       end if

       print '(A)', head
    else
       print '(A)', trim(line)
    end if

    if (guiding_centres) then
       do i=1,num_wann
          print '(2X, I0, ":s")', centers(i)
       end do
    end if

    if (a) then
       print '(A)', tail
    else
       call skip_to_end()
    end if
  end subroutine print_proj

  subroutine print_atoms(append, cart)
    logical, intent(in), optional :: append, cart
    logical                       :: a, c
    character(len=*), parameter   :: head = 'begin atoms_cart'
    character(len=*), parameter   :: tail = 'end atoms_cart'
    integer                       :: i

    a = .false.; c=.true.
    if (present(append)) a = append
    if (present(cart  )) c = cart

    appending: if (a) then
       print*
       print '(A)', head
       print '(A)', '  Bohr'
    else
       cartesian: if (c) then
          print '(A)', trim(line)
       else
          print '(A)', head
       end if cartesian

       read(*, '(A)') line
       unit: if (index(adjustl(lowercase(line)), 'bohr') == 1) then
          print '(A)', trim(line)
       else
          print '(A)', '  Bohr'
       end if unit
    end if appending

    do i=1,stru%nat
       print '(2X, I0, 3F9.5)', i, stru%pos(:, i) * stru%a
    end do

    if (a) then
       print '(A)', tail
    else
       call skip_to_end(printend=c)

       if (.not. c) print '(A)', tail
    end if
  end subroutine print_atoms

  subroutine print_kmesh(append)
    logical, intent(in), optional :: append
    logical                       :: a
    character(len=*), parameter   :: head = 'begin kpoints'
    character(len=*), parameter   :: tail = 'end kpoints'

    a = .false.
    if (present(append)) a = append

    if (a) then
       print*
       print '("!!! K-Points !!!")'
       print '("mp_grid : ", 3(I0, 1X))', mp_grid
       print '(A)', head
    else
       print '("mp_grid : ", 3(I0, 1X))', mp_grid
       print '(A)', trim(line)
    end if

    do i = 1, size(kmesh,1)
       print "(3F13.9)", kmesh(i, :)
    end do

    if (a) then
       print '(A)', tail
    else
       call skip_to_end()
    end if
  end subroutine print_kmesh

  subroutine print_kpath(append)
    logical, intent(in), optional :: append
    logical                       :: a
    character(len=*), parameter   :: head = 'begin kpoint_path'
    character(len=*), parameter   :: tail = 'end kpoint_path'

    a = .false.
    if (present(append)) a = append

    if (a) then
       print*
       print '("!!! BZ-Path for band structure !!!")'
       if (write_kpath) then
          i = findkey('bands_plot') ! must not print it twice
          if (i==0) print '(A)', 'bands_plot = .true.'
       end if
       print '(A)', head
    elseif (write_kpath) then    ! otherwise, print_block() will 
       print '(A)', trim(line)   ! include the ‘begin’
    end if

    if (write_kpath) then
       do i= 1, size(kpath,1)-1
          print "(A3,3F6.2,A4,3F6.2)", &
               trim(knames(i)),  kpath(i,:), &
               trim(knames(i+1)),kpath(i+1,:)
       enddo
    end if

    if (a) then
       print '(A)', tail
    elseif (write_kpath) then
       call skip_to_end()
    else
       call print_block()
    end if
  end subroutine print_kpath


!!!----------- Skip a generic block ---------
  subroutine skip_to_end(printend)
    logical, intent(in), optional :: printend
    logical                       :: p

    p = .true.
    if (present(printend)) p = printend

    skip: do
       if (index(adjustl(lowercase(line)), 'end') == 1) then
          if (p) print '(A)', trim(line)
          exit skip
       end if
       read(*, '(A)') line
    end do skip
  end subroutine skip_to_end


!!!----------- Reproduce a generic block ---------
  subroutine print_block()
    print '(A)', trim(line)

    block: do
       read(*, '(A)') line
       print '(A)', trim(line)
       if (index(adjustl(lowercase(line)), 'end') == 1) then
          exit block
       end if
    end do block
  end subroutine print_block


!!!----------- Key/value array management --------
  subroutine printkey(i)
    integer, intent(in) :: i

    print '(A, " = ", A)', trim(keys(i)), trim(vals(i))
  end subroutine printkey

  pure integer function findkey(key)
    character(len=*), intent(in) :: key

    do findkey = 1, nkeys
       if (lowercase(key) == keys(findkey)) goto 102
    end do

    findkey = 0

102 continue
  end function findkey

  subroutine allocate_keyval(sz)
    integer, intent(in) :: sz
    allocate(keys(sz), vals(sz))
  end subroutine allocate_keyval

  subroutine realloc_keyval(inc)
    integer, intent(in) :: inc
    character(len=BUFSZ), allocatable :: tmp(:)

    allocate(tmp(size(keys)+inc))
    tmp(1:size(keys)) = keys(:)
    call move_alloc(tmp, keys)

    allocate(tmp(size(vals)+inc))
    tmp(1:size(vals)) = vals(:)
    call move_alloc(tmp, vals)
  end subroutine realloc_keyval
end program write_win

!!/---
!! Local Variables:
!! mode: f90
!! End:
!!\---
