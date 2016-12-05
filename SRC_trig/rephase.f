!!! wien2wannier/SRC_trig/rephase.f
!!!
!!!    Rewrites inwf input file such that the phases given by psiargs
!!!    files are considered and the resulting Wannier functions should
!!!    contain less complex phases algorithm: the most probable actual
!!!    phase of a Wannier function is determined by parting the
!!!    interval [0,pi] in 100 slices and counting in which slice the
!!!    phases read in from the psiarg fall into. The slice with the
!!!    maximal count is then declared the "mean" phase. Note that
!!!    small phases are neglected due to possible numerical errors.
!!!
!!! Copyright 2009-2012 Philipp Wissgott
!!!           2013-2015 Elias Assmann

PROGRAM rephase
  use util,  only: line_count
  use const, only: DPk, TAU

  implicit none
  character(50) seedname


  logical :: flagfirst,dowrite
  integer iarg,argcount !command line input argument counter
  integer :: nemin,nemax
  integer :: Ljmax  !plain wave expansion cutoff
  integer :: nproj  !#projection = # final Wannier functions
  integer, allocatable :: norbitals(:)
  integer, allocatable :: lindex(:,:),mindex(:,:)      !l,m quantum number for orbitals
  integer :: j1,j2,j3,j            !counter
  integer,allocatable :: atomidx(:,:)
  integer :: widx ! index of wannier function to be rephased
  !character*1 orbital_character
  real(DPk), allocatable :: coeff(:,:,:)  !entry coefficients
  complex(DPk) :: ccoeff(5) !entry coefficients, complex version
  character(70) dummy,startmessage
  character(8) inwffileend
  character(9) psiargfileend
  real(DPk), allocatable :: phases(:,:)!phases
  integer, allocatable :: iphases(:) !integer index version of phases int(phases*100)
  integer :: maxel ! element in phase mesh with maximal occurance
  real(DPk) :: mphase !median phase from psiarg file
  integer, allocatable :: pcount(:)!counter of positive elements of phases
  integer :: nlines=0,pcounttmp
  character(2) :: idx
  integer :: counter(314) !counts values of the positive phase in int(pi*100)

  !default fileending: non spin-polarized
  inwffileend = ".inwf"
  psiargfileend = ".psiarg"
  startmessage = "++ Rewriting a non spin-polarized input file ++"

  !command line argument read-in
  iarg=command_argument_count()
  argcount = 1
  widx = 0
  dowrite = .false.
  if(iarg.ge.1) then
     do j=1,iarg
        call get_command_argument(j,dummy)
        if (dummy(1:1).eq.'-') then
           if ((dummy(2:3).eq.'up').or.(dummy(2:3).eq.'dn')) then
              !for spin-polarized calc. the fileendings have additional up/dn
              inwffileend = ".inwf"//dummy(2:3)
              psiargfileend = ".psiarg"//dummy(2:3)
              startmessage = "++ Rewriting a spin-polarized input file:"//dummy(2:3)//" ++"
           elseif (dummy(2:4).eq.'wf=') then
              !if this option is given, only the wannier function following the flag is
              !rephased, otherwise all Wannier function are rephased
              read(dummy(5:6),"(I2)")widx
           elseif (dummy(2:2).eq.'w') then
              !if this option is given inwf is rewritten,
              !otherwise the mean phase is only printed out
              dowrite = .true.
           elseif (dummy(2:2).eq.'h') then
              write(*,*)'finds the mean phases of Wannier functions reading'
              write(*,*)'the psiarg files after plotting.'
              write(*,*)'Usage: rephase [-up/-dn] [-w] [-wf m] case'
              write(*,*)'if the option "w" is given the inwf file is rewritten'
              write(*,*)'to take the phase into account(makeing any future'
              write(*,*)'Wannier function computed with that input file more'
              write(*,*)'likely to be real.'
              write(*,*)'with the "wf m" option only the Wannier function m is taken'
              write(*,*)'into account.'
              write(*,*)'Algorithm: the most probable actual phase of a Wannier function is'
              write(*,*)'determined by parting the interval [0,pi] in 100 slices'
              write(*,*)'counting in which slice the phases read in from the psiarg'
              write(*,*)'fall into. The slice with the maximal count is then declared'
              write(*,*)'the "mean" phase. Note that small phases are neglected due to '
              write(*,*)'possible numerical errors'
              stop
           else
              write(*,*)"Error: unknown option"
              stop
           endif
        else
           if (argcount.eq.1) then
              read(dummy,*)seedname
              argcount = argcount + 1
           else
              write(*,*)"Error: unknown option/input"
              stop
           endif
        endif
     enddo

  else
     write(*,*) 'Usage: rephase <seedname>'
     write(*,*) ' Options: -nw do not rewrite inwf, just compute phase(s)'
     write(*,*) ' Options: -w<n> compute phase only for nth Wannier function'
     stop
  endif

  if (dowrite) then
     write(*,*)startmessage
  endif

  open(unit=11,file=trim(seedname)//trim(inwffileend),status='old',action='read')

  read(11,"(A4)")dummy
  read (11,"(2I4,A60)") nemin,nemax,dummy
  read(11,"(2I3,A60)")Ljmax,nproj,dummy
  !allocate arrays
  allocate(norbitals(nproj))
  allocate(atomidx(nproj,5),lindex(nproj,5),mindex(nproj,5),coeff(nproj,5,2))
  allocate(pcount(nproj))

  write(*,*)"read-in initial values for Wannier functions..."
  flagfirst = .true.!first wannier function input to be changed
  do j1=1,nproj
     if ((widx.eq.j1).or.(widx.eq.0)) then
        !read-in psiarg file to compute median phase
        write(idx,"(I2)")j1
        open(unit=15,file=trim(seedname)//'_'//trim(idx) //trim(psiargfileend),status='old')
        if (flagfirst) then
           nlines = line_count(15)
           allocate(phases(nproj,nlines*10))
           flagfirst = .false.
        endif
        pcount(j1) = 0
        do j2=1,nlines
           read(15,"(10E16.8)") phases(j1,(j2-1)*10+1:(j2-1)*10+10)
           do j3=1,10
              !count positive entries
              if (phases(j1,(j2-1)*10+j3).gt.2.5d-1) then
                 pcount(j1) = pcount(j1) + 1
              endif
           enddo
        enddo
        close(15)
     endif

     read(11,"(I1,A60)")norbitals(j1),dummy
     do j2=1,norbitals(j1)
        read(11,"(2I2,I3,2F13.8,A47)")atomidx(j1,j2),lindex(j1,j2),mindex(j1,j2),coeff(j1,j2,1),coeff(j1,j2,2),dummy
     enddo
  enddo
  close(11)

  !reopen to write
  open(unit=11,file=trim(seedname)//trim(inwffileend),status='old',action='write')
  if (dowrite) then
     write(11,"(A4)")"BOTH"
     write(11,"(2I4,A60)")nemin,nemax,      "     # min band, max band                                 "
     write(11,"(2I3,A60)")Ljmax,nproj,"     # LJMAX max in exp(ibr) expansion, #Wannier functions"
  endif
  do j1=1,nproj
     ccoeff(:) = cmplx(0d0,0d0, DPk)
     if (dowrite) then
        write(11,"(I1,A60)")norbitals(j1),"     #number of entries per orbital (list of projected orbitals)"
     endif
     if ((widx.eq.j1).or.(widx.eq.0)) then
        allocate(iphases(pcount(j1)))
        pcounttmp = 1
        do j2=1,nlines
           do j3=1,10
              !set positive entries
              if (phases(j1,(j2-1)*10+j3).gt.2.5d-1) then
                 !integer conversion to index of phase mesh
                 iphases(pcounttmp) = mod(int(phases(j1,(j2-1)*10+j3)*50),int(TAU*50)) + 1
                 pcounttmp = pcounttmp + 1
              endif
           enddo
        enddo
        counter(:) = 0
        do j2=1,pcounttmp-1
           !count phases falling in certain slices of the phase mesh
           counter(iphases(j2)) = counter(iphases(j2)) + 1
        enddo
        !find the slice with maximal occurance
        maxel = 1
        do j2=2,314
           if (counter(j2).gt.counter(maxel)) then
              maxel = j2
           endif
        enddo
        !return to real values
        mphase = real(maxel)/100d0
        !discard small phases
        if ((mphase.ge.3.13d0).or.(mphase.le.0.02d0)) then
           mphase = 0d0
        endif
        !correct phase
        do j2=1,norbitals(j1)
           ccoeff(j2) = cmplx(coeff(j1,j2,1), coeff(j1,j2,2), DPk) * &
                & exp(-(0._DPk,1._DPk)*mphase)
        enddo
        write(*,*)"Wannier function:",j1," median phase:",mphase
     else
        !unchanged phase
        do j2=1,norbitals(j1)
           ccoeff(j2) = cmplx(coeff(j1,j2,1),coeff(j1,j2,2), DPk)
        enddo
        !write(*,*)"write init. coeff., WF:",j1," unchanged"

     endif
     if (dowrite) then
        do j2=1,norbitals(j1)
           write(11,"(2I2,I3,2F13.8,A47)")atomidx(j1,j2),lindex(j1,j2),mindex(j1,j2),&
                & real(ccoeff(j2)),aimag(ccoeff(j2))," # index of atom, L, M, coefficient (complex)"
        enddo
     endif

     if ((widx.eq.j1).or.(widx.eq.0)) then
        deallocate(iphases)
     endif

  enddo
  if (dowrite) then
     close(11)
  endif
  deallocate(norbitals)
  deallocate(atomidx,lindex,mindex,coeff)
  deallocate(pcount)

END PROGRAM rephase


!!/---
!! Local Variables:
!! mode: f90
!! End:
!!\---
!!
