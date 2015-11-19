PROGRAM find_bands
 !program finds the correct bands which are within an energy interval
 
 use util 

 implicit none

 integer :: unitout1,unitscf2 !file units
 character*50 seedname

 integer nkpoints,j1,j2
 integer iarg !cmmand line input argument counter
 real*4 Emin,Emax,tmp
 real*8 Fermienergy
 real*8 Eigenvalues(5)
 real*8 Eigenvalues_tot(10000)
 logical :: Ffound_Fermienergy
 integer nlines   !# lines in a file
 character*70 dummy,argdummy,startmessage
 integer :: counter,kidx,j
 logical ::  kread
 integer eigcount,eigcount_tot
 character*10 form
 character*10 output1fileend
 character*7 scf2fileend
 integer :: argcount

 !default fileending: non spin-polarized
 output1fileend = ".output1"
 scf2fileend  = ".scf2"
 startmessage = "++ find energy interval with standard files ++"
 
 !command line argument read-in
 iarg=iargc()
 argcount = 1
 if(iarg.ge.1) then
    do j=1,iarg
        call getarg(j,argdummy)
        if ((argdummy(1:1).eq.'-').and.(argcount.eq.1)) then
            if ((argdummy(2:3).eq.'up').or.(argdummy(2:3).eq.'dn')) then     
               !for spin-polarized calc. the fileendings have additional up/dn
               output1fileend = ".output1"//argdummy(2:3)
               scf2fileend  = ".scf2"//argdummy(2:3)               
               startmessage = "++ find energy interval with spin-polarized input files:"//argdummy(2:3)//" ++"
            elseif ((argdummy(2:5).eq.'soup').or.(argdummy(2:5).eq.'sodn')) then  
               output1fileend = ".outputso"
               scf2fileend  = ".scf2"//argdummy(4:5)               
               startmessage = "++ find energy interval with spin-polarized spin-orbit input files:"//argdummy(4:5)//" ++"
            elseif ((argdummy(2:2).eq.'h')) then
               write(*,*) 'finds the correct bands which are within an energy interval'
               write(*,*) 'Usage: find_bands [-up/-dn/-soup/-sodn] case Emin Emax (in eV)'
               stop
            else
               write(*,*) 'Unknown option'
               write(*,*) 'Usage: find_bands [-up/-dn/-soup/-sodn] case Emin Emax (in eV)'
               stop 
            endif
         else
            if (argcount.eq.1) then
               read(argdummy,*)seedname
               argcount = argcount + 1
               call getarg(j+1,dummy)
               read(dummy,*)Emin
               call getarg(j+2,dummy)
               read(dummy,*)Emax
            else
               cycle
            endif
        endif
     enddo
 else
   write(*,*) 'Usage: find_bands [-up/-dn/-soup/-sodn] case Emin Emax (in eV)'
   stop 
 endif

  write(*,*)startmessage

  unitout1 = 12
  unitscf2 = 11

  !get Fermi energy
  open(unit=unitscf2,file=clearspace(seedname)//clearspace(scf2fileend),status='old')
  nlines=line_count(unitscf2)
  Ffound_Fermienergy = .false.
  do j1=1,nlines
     read(unitscf2,"(A60)")dummy
     if (dummy(1:4).eq.":FER") then
       Ffound_Fermienergy = .true.
       read(dummy(39:48),"(F10.5)")Fermienergy
     endif
  enddo
  close(unitscf2)
  if (Ffound_Fermienergy.eqv. .false.) then
    write(*,*)"Did not find Fermi energy in "//clearspace(seedname)//clearspace(scf2fileend)
    stop
  else
     write(*,*)"Fermi energy(from "//clearspace(seedname)//clearspace(scf2fileend)//&
               &" file) is: ",sngl(Fermienergy),"Ry=",sngl(Fermienergy*13.6),"eV"
  endif
  
  !open seed.output1
  open(unit=unitout1,file=clearspace(seedname)//clearspace(output1fileend),status='old')
  nlines=line_count(unitout1)
  counter = 0
  kidx = 0
  kread =.false.
  do while (counter.lt.nlines)
     counter = counter  + 1
     form = "(A70)"
     read(unitout1,form)dummy
     if ((dummy(15:42).eq."EIGENVALUES BELOW THE ENERGY")) then !.or.(clearspace(dummy).eq."")) then
        !turn the eigenvalue read off
        kread = .false.
        call print_eigenvalues_information(kidx,Eigenvalues_tot,eigcount_tot,Emin,Emax,Fermienergy)
     endif 
     if (kread.eqv..true.) then
        Eigenvalues(:) = 0d0
        !how many eigenvalues are there in this line
        eigcount = char_count(clearspace(dummy(8:8)//dummy(21:21)//dummy(34:34)//dummy(47:47)//dummy(60:60)))
        if (eigcount.eq.0) then
          cycle
        endif
        write(form,"(A1,I1,A6)")'(',eigcount,'F13.7)'
        read(dummy(5:70),form)Eigenvalues_tot(eigcount_tot+1:eigcount_tot+eigcount)!Eigenvalues(1:eigcount)
        eigcount_tot = eigcount_tot + eigcount
     endif
     
     
     if (dummy(6:20).eq."EIGENVALUES ARE") then
        !turn the eigenvalue read on
        kidx = kidx + 1
        eigcount_tot = 0
        Eigenvalues_tot(:) = 0d0
        kread = .true.
     endif
  enddo
  close(unitout1)
  
END PROGRAM

SUBROUTINE print_eigenvalues_information(kidx,Eigenvalues_tot,eigcount_tot,Emin,Emax,Fermienergy)

implicit none

 real*4 Emin,Emax
 real*8 Fermienergy
 real*8 Eigenvalues_tot(1000)
 integer :: kidx
 integer :: eigcount_tot
 integer :: j1
 integer minidx,maxidx
 logical flagmin
 
 !consider Fermi energy and convert to eV
 flagmin = .false.
 do j1=1,eigcount_tot
   Eigenvalues_tot(j1) = Eigenvalues_tot(j1) - Fermienergy
   Eigenvalues_tot(j1) = Eigenvalues_tot(j1)*13.6d0
   if ((Eigenvalues_tot(j1).ge.Emin).and.(flagmin.eqv..false.))then
      minidx = j1
      flagmin = .true.
   endif
  if (Eigenvalues_tot(j1).le.Emax) then
      maxidx = j1
   endif
 enddo
 write(*,"(A2,I4,A9,I3,A15,I3,A15,I3,A22,I3)")"k=",kidx," | #eig.=",eigcount_tot," | min in int.=",minidx,&
                                     " | max in int.=",maxidx," | #bands in en.-int.=",maxidx-minidx+1

END SUBROUTINE
