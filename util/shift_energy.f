PROGRAM shift_energy
 !program shifts the energy within the seed.eig file to account 
 !for the Fermi energy
 !P.Wissgott 10/01/10
 
 use util 

 implicit none

 integer :: uniteig,unitscf2 !file units
 character*50 seedname

 integer iarg !cmmand line input argument counter
 real*8 Fermienergy
 real*8, allocatable :: Energies(:)
 integer,allocatable :: bandidx(:),kidx(:)
 logical :: Ffound_Fermienergy,userinput
 integer nlines   !# lines in a file
 integer :: j1,j !loop variable
 character*70 dummy,startmessage
 character*70 argdummy
 integer :: argcount
 character*6 eigfileend
 character*7 scf2fileend

  !default fileending: non spin-polarized
  scf2fileend  = ".scf2"
  eigfileend = ".eig"

  startmessage = "++ Shifting energy in standard eig file ++"
  userinput = .false.

 !command line argument read-in
 iarg=iargc()
 argcount = 1
 if(iarg.ge.1) then
     do j=1,iarg
        call getarg(j,argdummy)
        if (argdummy(1:1).eq.'-') then
           if ((argdummy(2:3).eq.'up').or.(argdummy(2:3).eq.'dn')) then     
              !for spin-polarized calc. the fileendings have additional up/dn
              scf2fileend  = ".scf2"//argdummy(2:3)
              eigfileend = ".eig"//argdummy(2:3)
              startmessage = "++ Shifting energy in eig file: "//argdummy(2:3)//" ++"
           elseif (argdummy(2:3).eq.'ef') then     
               userinput = .true.
               call getarg(j+1,argdummy)
               read(argdummy,"(F5.3)")Fermienergy 
               Fermienergy = Fermienergy*13.6
               userinput = .true.
           elseif (argdummy(2:2).eq.'h') then
               write(*,*) 'program shifts the energy within the case.eig file to account for the Fermi energy'
               write(*,*) 'read from the scf2[up/dn] file'
               write(*,*) 'Usage: shift_energy [-up/-dn] [-ef EF] case'
               write(*,*) 'if the option "-ef EF" is given, the shifting energy is set to EF(in eV)'
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

 uniteig = 12
 unitscf2 = 11

  !get Fermi energy
  if (.not. userinput) then
     open(unit=unitscf2,file=clearspace(seedname)//clearspace(scf2fileend),status='old')
     nlines=line_count(unitscf2)
     Ffound_Fermienergy = .false.
     do j1=1,nlines
       read(unitscf2,"(A60)")dummy
       if (dummy(1:4).eq.":FER") then
         Ffound_Fermienergy = .true.
         !write(*,*)dummy(39:48)
         read(dummy(39:48),"(F10.5)")Fermienergy
         !write(*,*)dummy
       endif
     enddo
     close(unitscf2)
     if (Ffound_Fermienergy.eqv. .false.) then
       write(*,*)"Did not find Fermi energy in "//clearspace(seedname)//clearspace(scf2fileend)
       stop
     endif
     !conversion to eV
     Fermienergy = Fermienergy*13.6
     write(*,*)"Fermi energy(from "//clearspace(seedname)//clearspace(scf2fileend)//" file) is: ", & 
                     &sngl(Fermienergy/13.6),"Ry       (=",sngl(Fermienergy),"eV)"
  else
     write(*,*)"Fermi energy(from user input) is: ",sngl(Fermienergy/13.6),"Ry     (=",sngl(Fermienergy),"eV)"
  endif
  

  !open seed.eig
  open(unit=uniteig,file=clearspace(seedname)//clearspace(eigfileend),status='old')
  nlines=line_count(uniteig)
  allocate(Energies(nlines))
  allocate(bandidx(nlines),kidx(nlines))
  do j1=1,nlines
      read(uniteig,"(2I12,F22.16)")bandidx(j1),kidx(j1),Energies(j1)
  enddo
  close(uniteig)
  
  !reopen seed.eig to write
  open(unit=uniteig,file=clearspace(seedname)//clearspace(eigfileend),status='old')
  do j1=1,nlines
      write(uniteig,"(2I12,F22.16)")bandidx(j1),kidx(j1),Energies(j1)-Fermienergy
  enddo
  close(uniteig)

  deallocate(Energies,bandidx,kidx)
  
END PROGRAM