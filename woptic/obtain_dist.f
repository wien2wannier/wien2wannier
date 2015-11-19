program obtain_dist

implicit none

character*80 seed,seed2,arg,dummy
integer :: unit1,unit2,nk,nw,j1,j2,j3,co,stat
real*8,allocatable :: centres(:,:),distmatx(:,:),distmaty(:,:),distmatz(:,:)
real*8,allocatable :: r(:,:)
real*8 :: normus1,normus2,U,dist,tmp(10),factor

  unit1 = 1
  unit2 = 2

  call get_command_argument(1, VALUE=arg, STATUS=stat)
  if (stat /= 0) stop 'error getting seed'
  read(arg, *) seed
  write(*,*)seed
  
  !obtain number of Wanniers
!   open(UNIT=unit1, FILE=trim(seed)//".ham", STATUS='old', FORM='formatted')
!   read(unit1,*)nk,nw
!   write(*,*)nk,nw
!   close(unit1)
  open(UNIT=unit1, FILE=trim(seed)//"_centres.xyz", STATUS='old', FORM='formatted')
  read(unit1,*)nw
  read(unit1,*)
  allocate(centres(nw,3),distmatx(nw,nw),distmaty(nw,nw),distmatz(nw,nw))
  do j1=1,nw
     read(unit1,*)dummy,centres(j1,:)
     write(*,*)centres(j1,:)
  enddo
  close(unit1)

  do j1=1,nw
     do j2=1,nw
        distmatx(j1,j2) = 1.889725d0*(centres(j1,1)-centres(j2,1))
        distmaty(j1,j2) = 1.889725d0*(centres(j1,2)-centres(j2,2))
        distmatz(j1,j2) = 1.889725d0*(centres(j1,3)-centres(j2,3))
     enddo
  enddo

  open(UNIT=unit2, FILE=trim(seed)//".distmatrix", STATUS='unknown', FORM='formatted')
  do j1=1,nw
     write(unit2,"(100F12.4)")distmatx(j1,:)
  enddo
  do j1=1,nw
     write(unit2,"(100F12.4)")distmaty(j1,:)
  enddo
  do j1=1,nw
     write(unit2,"(100F12.4)")distmatz(j1,:)
  enddo
  close(unit2)

  

end program