!!! wien2wannier/SRC_trig/vec2ascii.f
!!!
!!!    Translates WIEN2k vector files to plain text.  Based on
!!!    join_vectorfiles.
!!!
!!!    Usage: vec2ascii [-up/-dn] [-c] <case> <numberofparallelfiles>
!!!
!!! Copyright 2013-2015 Elias Assmann

PROGRAM vec2ascii
  use util,      only: line_count
  use structmod, only: struct_t, struct_read
  use const,     only: BUFSZ, R8, C16

  implicit none

 character(50) seedname
 integer nkpoints,nfiles
 integer iarg,argcount !command line input argument counter
 integer jatom,i,j,k,jl,jj,jk
 integer lmax,lomax,nloat
 integer unitklist,unitkgen,unitstruct,unitvector, unitenergy, unittargetvector,unittargetenergy
 character(70) argdummy,startmessage
 logical :: usecomplex
 INTEGER            NE, NV
 DOUBLE PRECISION   SX, SY, SZ, WEIGHT
 CHARACTER(3)       IPGR,filenr
 CHARACTER(10)      KNAME
 DOUBLE PRECISION, allocatable ::  EIGVAL(:)
 DOUBLE PRECISION, allocatable ::   E(:,:),ELO(:,:,:)
 INTEGER, pointer :: KZZ(:,:)
 REAL(R8), allocatable ::  Z(:,:)
 COMPLEX(C16), allocatable ::  ZC(:,:)
 type(struct_t) stru
 character(bufsz) :: vectorfileend, targetvectorfileend
 character(bufsz) :: energyfileend, targetenergyfileend
 character(bufsz) :: suf
 DOUBLE PRECISION   eorb_ind


  !default fileending: non spin-polarized
  targetvectorfileend = ".vector_ascii"
  vectorfileend = ".vector"
  energyfileend = ".energy"
  targetenergyfileend = ".energy_ascii"
  startmessage = "++ join vector files of standard input ++"
  unitkgen = 1
  unitvector = 2
  unitenergy = 3
  unittargetvector = 4
  unittargetenergy = 5
  unitklist = 11

  !parameters
  LMAX = 13
  LOMAX = 3
  nloat = 3

  !command line argument read-in
  iarg=command_argument_count()
  argcount = 1


  usecomplex = .false.
  if ((iarg.ge.1)) then
    do j=1,iarg
        call get_command_argument(j,argdummy)
        if (argdummy(1:1).eq.'-') then
            if ((argdummy(2:3).eq.'up').or.(argdummy(2:3).eq.'dn')) then
               !for spin-polarized calc. the fileendings have additional up/dn
              vectorfileend = ".vector"//argdummy(2:3)
              energyfileend = ".energy"//argdummy(2:3)
              startmessage = "++ join vector files of spin-polarized input files:"//argdummy(2:3)//" ++"
            elseif ((argdummy(2:5).eq.'soup').or.(argdummy(2:5).eq.'sodn')) then
               vectorfileend = ".vectorso"//argdummy(4:5)
               energyfileend = ".energyso"//argdummy(4:5)
               usecomplex = .true.
               startmessage = "++ ascii'ize vector files of spin-polarized spin-orbit input files:"//argdummy(4:5)//" ++"
            elseif (argdummy(2:2).eq.'c') then
               usecomplex = .true.
            elseif (argdummy(2:2).eq.'h') then
               write(*,*) 'joins multiple WIEN2K vector files to one for further processing'
               write(*,*) 'Usage: vec2ascii [-up/-dn/-soup/-sodn/-c] case numberofparallelfiles'
               stop
            else
               write(*,*) 'Unknown option'
               write(*,*) 'Usage: vec2ascii [-up/-dn/-soup/-sodn/-c] case numberofparallelfiles'
               stop
            endif
         else
            if (argcount.eq.1) then
               read(argdummy,*)seedname
               argcount = argcount + 1
            else
               read(argdummy,*)nfiles
            endif
        endif
     enddo
 else
   write(*,*) 'Usage: vec2ascii [-up/-dn/-soup/-sodn/-c] case numberofparallelfiles'
   stop
 endif

 write(*,*)"case is:", seedname
 write(*,*)startmessage

 open(unit=unitklist,file=trim(seedname)//'.klist',status='old')
 open(unit=unitstruct,file=trim(seedname)//'.struct',status='old')
 call struct_read(unitstruct, stru)
 nkpoints = line_count(unitklist) -2 !take care
 write(*,*)"#kpoints",nkpoints
 close(unitstruct)
 close(unitklist)


 allocate( E(LMAX,stru%nneq) )
 allocate( ELO(0:LOMAX,nloat,stru%nneq) )

 open(unit=unittargetvector,file=trim(seedname)//trim(targetvectorfileend), &
      & status='unknown',form='formatted')
 open(unit=unittargetenergy,file=trim(seedname)//trim(targetenergyfileend), &
      & status='unknown',form='formatted')

 do j=1,nfiles
    if (nfiles == 1) then
       suf = ""
    else
       write(filenr,"(I3)")j
       suf = "_"//trim(adjustl(filenr))
    end if

    open(unit=unitvector, &
         file=trim(seedname)//trim(vectorfileend)//suf, &
         status='old',form='unformatted')
    open(unit=unitenergy, &
         file=trim(seedname)//trim(energyfileend)//suf, &
         status='old',form='formatted')
    do jatom= 1,stru%nneq
       read(unitvector) (E(jl,jatom),jl=1,LMAX)
       read(unitvector) ((ELO(jl,k,jatom),jl=0,LOMAX),k=1,nloat)
       read(unitenergy,'(100(f9.5))') (E(jl,jatom),jl=1,LMAX),eorb_ind
       read(unitenergy,'(100(f9.5))') ((ELO(jl,k,jatom),jl=0,LOMAX),k=1,nloat)
       if (j.eq.1) then
          write(unittargetvector,*) (E(jl,jatom),jl=1,LMAX)
          write(unittargetvector,*) ((ELO(jl,k,jatom),jl=0,LOMAX),k=1,nloat)
          write(unittargetenergy,'(100(f9.5))') (E(jl,jatom),jl=1,LMAX),eorb_ind
          write(unittargetenergy,'(100(f9.5))') ((ELO(jl,k,jatom),jl=0,LOMAX),k=1,nloat)
       endif
    enddo
    do jk=1,nkpoints

       read(unitvector,end=888,err=888) SX, SY, SZ, KNAME, NV, NE, WEIGHT!, IPGR
       allocate( KZZ(3,NV) )
       read(unitenergy,'(3e19.12,a10,2i6,f5.1,a3)') SX, SY, SZ, KNAME, NV, NE, WEIGHT, IPGR
       read(unitvector) (KZZ(1,I),KZZ(2,I),KZZ(3,I),I=1,NV)
       write(unittargetvector,*) SX, SY, SZ, KNAME, NV, NE, WEIGHT!, IPGR
       write(unittargetenergy,'(3e19.12,a10,2i6,f5.1,a3)') &
            & SX, SY, SZ, KNAME, NV, NE, WEIGHT, IPGR
       write(unittargetvector,*) (KZZ(1,I),KZZ(2,I),KZZ(3,I),I=1,NV)
       allocate(Z(NV,NE),ZC(NV,NE))
       allocate(EIGVAL(NE))
       do jj = 1, NE
         read(unitvector) I, EIGVAL(I)
         write(unittargetvector,*) I,EIGVAL(I)
         read(unitenergy,*) I, EIGVAL(I)
         write(unittargetenergy,*) I,EIGVAL(I)
         if (usecomplex) then
             read(unitvector) (ZC(jl,I),jl=1,NV)
             write(unittargetvector,*) (ZC(jl,I),jl=1,NV)
         else
             read(unitvector) (Z(jl,I),jl=1,NV)
             write(unittargetvector,*) (Z(jl,I),jl=1,NV)
         endif
       enddo
      deallocate(Z,ZC,EIGVAL,KZZ)
 888  continue
    enddo
    close(unitvector)
    close(unitenergy)
 print*, 'reading energy/vector file',j,' done'
 enddo
 close(unittargetvector)
 close(unittargetenergy)
 write(*,*)"++ joining vector files(+energy files) complete ++"
END PROGRAM vec2ascii


!!/---
!! Local Variables:
!! mode: f90
!! End:
!!\---
