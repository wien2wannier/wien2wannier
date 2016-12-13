!!! wien2wannier/SRC_trig/write_win_backend.f
!!!
!!!    Prepares input case.win for Wannier90
!!!
!!! Copyright 2009-2012 Philipp Wissgott
!!!           2013-2015 Elias Assmann

program write_win
  use structmod, only: struct_t, struct_read
  use inwfmod,   only: inwf_t, inwf_read
  use const,     only: BUFSZ, DPk
  use util,      only: lowercase, newunit
  use clio,      only: fetcharg, argstr
  use kpoints,   only: get_kmesh_band, get_kmesh_klist
  use reallocate,only: realloc

  implicit none

!!!----------- Configurable parameters -----------
  !! Initial size and size increment of keys, vals.
  integer, parameter :: NKEY_FIRST=50, NKEY_INC=50

  character(len=*), parameter ::                       &
       fmt_kpoint_path  = "(2(A3,3(1X,F9.5),2X))",     &
       fmt_brlat        = '(3(1X,F11.6))',             &
       fmt_centers      = '(I4, ":s")',                &
       fmt_atoms        = '(I4, 3(1X,F13.8))',         &
       fmt_mp_grid      = '("mp_grid : ", 3(I0, 1X))', &
       fmt_mp_grid_bare = '(              3(I0, 1X))', &
       fmt_kmesh        = "(3F19.15)"

!!!----------- Variables               -----------
  integer :: num_bands, mp_grid(3)
  integer :: c, i, j, k, nfill, iarg

  real(DPk), allocatable :: kpath(:,:), kmesh(:,:)
  character, allocatable :: knames(:)
  integer,   allocatable :: centers(:)

  type(argstr) :: inwffile, structfile, klistfile, bandfile

  character(len=BUFSZ) :: line, key, filler, comment

  character(len=BUFSZ), allocatable :: keys(:), vals(:)
  logical,              allocatable :: keys_done(:)
  integer                           :: nkeys

  logical :: atoms_done=.false., uc_done=.false., proj_done=.false.
  logical :: kmesh_done=.false., kpath_done=.false., mpgrid_done=.false.
  logical :: bandsplot_done=.false., guiding_done=.false.
  logical :: write_kpath, bands_plot, guiding_centres

  type(struct_t) :: stru
  type(inwf_t)   :: inwf

!!!----------- Code                    -----------
  call fetcharg(1, inwffile,   "failed to get `inwf' argument")
  call fetcharg(2, structfile, "failed to get `struct' argument")
  call fetcharg(3, klistfile,  "failed to get `klist' argument")
  call fetcharg(4, bandfile,   "failed to get `klist_band' argument")

!!! Read ‘inwf’ file for num_bands, Nproj
  call inwf_read(inwffile, inwf)

  num_bands = inwf%bmax - inwf%bmin + 1

!!! Read parameters from command line
  call allocate_keyval(NKEY_FIRST)

  !! num_wann, num_bands, and mp_grid must match with ‘inwf’ and
  !! ‘klist’
        keys(1)        = 'num_bands'
  write(vals(1), '(I0)')  num_bands
  ! rationale for setting num_wann: num_wann=0 is not very useful, and
  ! this matches the behavior of w2w
        keys(2)        = 'num_wann'
  write(vals(2), '(I0)')  merge(inwf%Nproj, num_bands, inwf%Nproj>0)
  ! the following will be filled later
        keys(3) = 'mp_grid'
        vals(3) = ''

  nkeys = 3
  optarg: do iarg = 5, command_argument_count(), 2
     nkeys = nkeys+1
     if (nkeys > size(keys)) call realloc_keyval(NKEY_INC)

     call fetcharg(iarg,   keys(nkeys))
     i = findkey(keys(nkeys))
     if (i==nkeys) then
        i = nkeys
     else
        nkeys = nkeys-1
     end if
     call fetcharg(iarg+1, vals(i))
  end do optarg

!!! Find centers of projections for “guiding centres”; disable
!!! “guiding centres” if a projection is not uniquely centered, or if
!!! we have no projections
  allocate(centers(inwf%Nproj))

  guiding_centres = inwf%Nproj>0

  do i = 1, inwf%Nproj
     ylm: do j = 1, inwf%projections(i)%NY
        c = inwf%projections(i)%iat(j)
        if (j==1) then
           centers(i) = c
        elseif (c /= centers(i)) then
           guiding_centres = .false.
        end if
     end do ylm
  end do

  ! user input always overrides
  i = findkey('guiding_centres')
  if (i /= 0) read(vals(i), *) guiding_centres

!!! Read ‘struct’
  call struct_read(structfile%s, stru)

!!! Read ‘klist’ for MP k-mesh; read mp_grid from ‘klist’ unless it
!!! was given as an option
  i = findkey('mp_grid')
  if (vals(i)=='') then
     call get_kmesh_klist(klistfile%s, kmesh, stru, mp_grid)
     write(vals(i), fmt_mp_grid_bare) mp_grid
  else
     call get_kmesh_klist(klistfile%s, kmesh, stru)
  end if

!!! Read ‘klist_band’ for BZ path
  inquire(FILE=bandfile%s, EXIST=write_kpath)
  bands_plot = write_kpath
  i = findkey('bands_plot')
  if (i /= 0) read(vals(i), *) bands_plot
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

     special_keys: select case (key)
     case ('bands_plot')
        bandsplot_done = .true.
     case ('guiding_centres')
        guiding_done   = .true.
     case ('mp_grid')
        mpgrid_done    = .true.
     end select special_keys
  end do template
101 continue

!!! Add any parameters that were not in the template
  if (.not. all(keys_done)) print*

  rest: do i=1,nkeys
     if (keys_done(i)) cycle rest

     ! the following keys will be written together with the
     ! corresponding blocks if those have not yet appeared
     if (keys(i) == 'bands_plot'      .and. .not. kpath_done) cycle rest
     if (keys(i) == 'guiding_centres' .and. .not. proj_done)  cycle rest
     if (keys(i) == 'mp_grid'         .and. .not. kmesh_done) cycle rest

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
       print '(A)', 'Bohr'
    else
       print '(A)', trim(line)

       read(*, '(A)') line
       unit: if (index(adjustl(lowercase(line)), 'bohr') == 1) then
          print '(A)', trim(line)
       else
          print '(A)', 'Bohr'
       end if unit
    end if appending

    print fmt_brlat, stru%prim_dir

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
       if (.not. guiding_done) &
            print '("guiding_centres = ", L1)', guiding_centres
       print '(A)', head
    else
       print '(A)', trim(line)
    end if

    do i=1,inwf%Nproj
       print fmt_centers, centers(i)
    end do

    if (a) then
       print '(A)', tail
    else
       call skip_to_end()
    end if
  end subroutine print_proj

  subroutine print_atoms(append, cart)
    logical, intent(in), optional :: append, cart
    logical                       :: a, c
    character(len=*), parameter   :: &
         head_c = 'begin atoms_cart', tail_c = 'end atoms_cart', &
         head_f = 'begin atoms_frac', tail_f = 'end atoms_frac'
    integer                       :: i

    a = .false.; c=.true.
    if (present(append)) a = append
    if (present(cart  )) c = cart

    appending: if (a) then
       print*
       if (c) then
          print '(A)', head_c
          print '(A)', 'Bohr'
       else
          print '(A)', head_f
       end if
    else
       print '(A)', trim(line)

       if (c) then
          read(*, '(A)') line
          unit: if (index(adjustl(lowercase(line)), 'bohr') == 1) then
             print '(A)', trim(line)
          else
             print '(A)', 'Bohr'
          end if unit
       end if
    end if appending

    do i=1,stru%nat
       if (c) then
          print fmt_atoms, i, matmul(stru%conv_dir, stru%pos(:, i))
       else
          print fmt_atoms, i, matmul(transpose(stru%stru2frac), stru%pos(:, i))
       end if
    end do

    if (a) then
       if (c) then
          print '(A)', tail_c
       else
          print '(A)', tail_f
       end if
    else
       call skip_to_end()
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
       if (.not. mpgrid_done) print fmt_mp_grid, mp_grid
       print '(A)', head
    else
       print '(A)', trim(line)
    end if

    do i = 1, size(kmesh,1)
       print fmt_kmesh, kmesh(i, :)
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
       if (.not. bandsplot_done) print '("bands_plot = ", L1)', bands_plot
       print '(A)', head
    elseif (write_kpath) then    ! otherwise, print_block() will
       print '(A)', trim(line)   ! include the ‘begin’
    end if

    if (write_kpath) then
       do i= 1, size(kpath,1)-1
          print fmt_kpoint_path, trim(knames(i)),  kpath(i,:), &
               &                 trim(knames(i+1)),kpath(i+1,:)
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
       if (lowercase(key) == keys(findkey)) return
    end do

    findkey = 0
  end function findkey

  subroutine allocate_keyval(sz)
    integer, intent(in) :: sz
    allocate(keys(sz), vals(sz))
  end subroutine allocate_keyval

  subroutine realloc_keyval(inc)
    integer, intent(in) :: inc

    call realloc(keys, shape(keys)+inc)
    call realloc(vals, shape(vals)+inc)
  end subroutine realloc_keyval
end program write_win


!!/---
!! Local Variables:
!! mode: f90
!! End:
!!\---
