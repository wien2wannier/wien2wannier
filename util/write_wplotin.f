PROGRAM write_win
 !program prepares input seed.wplotin for wplot
 ! P.Wissgott 10/01/10
 
 use util 

 implicit none
 character*50 :: seedname
 logical :: writeout
 integer iarg !command line input argument counter
 integer :: origin(3),divisor
 integer vec(3)
 integer :: nlines   !# lines in a file
 integer :: j1,j,argcount
 character*80 dummy,argdummy,startmessage
 character*10 wplotinfileend
 character*7 woutfileend
 

  wplotinfileend = ".wplotin"
  woutfileend = ".wout"
 startmessage = "++ Preparing a standard input file ++"

  !command line argument read-in
 iarg=iargc()
 argcount = 1
 if(iarg.ge.1) then
     do j=1,iarg
        call getarg(j,argdummy)
        if (argdummy(1:1).eq.'-') then
           if ((argdummy(2:3).eq.'up').or.(argdummy(2:3).eq.'dn')) then     
              !for spin-polarized calc. the fileendings have additional up/dn
              wplotinfileend = ".wplotin"//argdummy(2:3)
              woutfileend = ".wout"//argdummy(2:3)
              startmessage = "++ Preparing a spin-polarized input file:"//argdummy(2:3)//" ++"
           elseif (argdummy(2:2).eq.'h') then
               write(*,*) 'prepares input file case.wplotin[up/dn] for wplot'
               write(*,*) 'Usage: write_wplotin [-up/-dn] case'
               stop 
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
  endif

  write(*,*)startmessage
  write(*,*)"From wannier90 output "//clearspace(seedname)//clearspace(woutfileend)&
                                    & //"(in Angstrom):"

  !read wannier90 output
  open(unit=1,file=clearspace(seedname)//clearspace(woutfileend),status='old')
  nlines = line_count(1)
  writeout = .false.
  do j1=1,nlines
      read(1,"(A80)")dummy
      if (dummy(2:12).eq."Final State") then
         writeout = .true.
      elseif (dummy(2:12).eq."           ") then
         writeout = .false.
      endif
      if (writeout.eqv..true.) then
        write(*,*)dummy
      endif
  enddo
  close(1)
  

  !write the input file
  open(unit=1,file=clearspace(seedname)//clearspace(wplotinfileend),status='unknown')
  write(1,"(A52)")"3D ORTHO        # mode O(RTHOGONAL)|N(ON-ORTHOGONAL)"
  print*,"Enter origin of spatial mesh in fractions of the conv. unit vectors,[n1 n2 n3 ndiv]:"
  read (*,*) origin(:),divisor
  write(1,"(I2,2I3,I2,A30)")origin(:),divisor,"     #x, y, z, divisor of origin"
  print*,"Enter endpoint of spatial mesh on x-axis in frac. of conv. u. vec.,[n1 n2 n3 ndiv]:"
  read (*,*) vec(:),divisor
  write(1,"(I2,2I3,I2,A31)")vec(:),divisor,"    #x, y, z, divisor of x-end"
  print*,"Enter endpoint of spatial mesh on y-axis in frac. of conv. u. vec.,[n1 n2 n3 ndiv]:"
  read (*,*) vec(:),divisor
  write(1,"(I2,2I3,I2,A31)")vec(:),divisor,"    #x, y, z, divisor of y-end"
  print*,"Enter endpoint of spatial mesh on z-axis in frac. of conv. u. vec.,[n1 n2 n3 ndiv]:"
  read (*,*) vec(:),divisor
  write(1,"(I2,2I3,I2,A31)")vec(:),divisor,"    #x, y, z, divisor of z-end"
  print*,"Enter number of mesh points,[nx ny nz]:"
  read (*,*) vec(:)
  write(1,"(3I3,A40)")vec(:)," 0 0 0 # grid points and echo increments"
  write(1,"(A50)")"NO              # DEP(HASING)|NO (POST-PROCESSING)"
  write(1,"(A48)")"WAN ANG LARGE   # switch ANG|ATU|AU LARGE|SMALL   "
  write(1,"(A40)")"1  1            # k-point, Wannier index             "

 ! read(1,*)dummy
 ! read(1,"(2I3,A60)")nemin,nemax,dummy
 ! options%num_bands = nemax - nemin + 1
 ! read(1,"(2I3,A60)")tmpint,options%num_wann
  close(1)

  !print*,"Enter minimal and maximal Bloch band,[n1 n2]:"

 

END PROGRAM



