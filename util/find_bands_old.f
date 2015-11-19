PROGRAM find_bands
 !program finds the correct bands which are within an energy interval
 
 use util 

 implicit none

 integer :: unitout1,unitscf2 !file units
 character*50 seedname

 integer line_count,nkpoints,j1,j2
 integer iarg,j,argcount !cmmand line input argument counter
 real*4 Emin,Emax,tmp
 real*8 Fermienergy
 real*8 Eigenvalues(5)
 real*8 Eigenvalues_tot(1000)
 logical :: Ffound_Fermienergy
 !integer jelement,jmult,jwrite !loop variables
 !integer mult !multiplicity
 !integer atomcount !counter for unit cell atoms
 integer nlines   !# lines in a file
 !integer nkpoints_band !#kpoints for bandstructure plots
 !real*8, allocatable :: kpoints(:,:),kpoints_band(:,:) !k-mesh points
 !character*6, allocatable :: kpoints_band_spec(:)
 !real*8 :: tmp2(3)
 character*70 dummy
 integer :: counter,kidx
 logical ::  kread
 integer eigcount,eigcount_tot
 character*10 form
 character*7 scf2fileend
 character*10 output1fileend

 !default fileending: non spin-polarized
 scf2fileend = ".scf2  "
 output1fileend = ".output1"

 iarg=iargc()
 argcount = 1
 if(iarg.ge.3) then
     do j=1,iarg
        call getarg(j,dummy)
        if (dummy(1:2).eq.'--') then
           if ((dummy(3:4).eq.'up').or.(dummy(3:4).eq.'dn')) then     
              !for spin-polarized calc. the fileendings have additional up/dn
              scf2fileend = ".scf2"//dummy(3:4)
              output1fileend = ".output1"//dummy(3:4)
              !write(*,*)"changing the end"
           else
              write(*,*)"Error: unknown option"
              stop
           endif
        else
            if (argcount.eq.1) then
               read(dummy,*)seedname
               argcount = argcount + 1
            elseif (argcount.eq.2) then
               read(dummy,*)Emin   
               argcount = argcount + 1
            elseif (argcount.eq.3) then
               read(dummy,*)Emax
               argcount = argcount + 1
            else
               write(*,*)"Error: unknown option/input"
               stop
            endif
        endif
     enddo
  
  else
     write(*,*) 'Error: At least 3 commandline arguments have to be given.'
     stop
  endif
  if (Emin.ge.Emax) then
     tmp = Emin
     Emin = Emax
     Emax = tmp
  endif
  write(*,*)"debug",seedname,Emin,Emax,scf2fileend

  unitout1 = 12
  unitscf2 = 11

  !get Fermi energy
  open(unit=unitscf2,file=clearspace(seedname)//clearspace(scf2fileend),status='old')
  nlines=line_count(clearspace(seedname)//clearspace(scf2fileend),unitscf2)
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
  else
     write(*,*)"Fermi energy(from "//clearspace(seedname)//clearspace(scf2fileend)//" file) is: ",&
                                                 &sngl(Fermienergy),"Ry=",sngl(Fermienergy*13.6),"eV"
  endif
  

  !open seed.output1
  open(unit=unitout1,file=clearspace(seedname)//clearspace(output1fileend),status='old')
  nlines=line_count(clearspace(seedname)//clearspace(output1fileend),unitout1)
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
        !write(*,*)"off"
        if (dummy(15:42).eq."EIGENVALUES BELOW THE ENERGY") then
          call print_eigenvalues_information(kidx,Eigenvalues_tot,eigcount_tot,Emin,Emax,Fermienergy)
        endif
     endif 
     if (kread.eqv..true.) then
        !write(*,*)dummy
        Eigenvalues(:) = 0d0
        !how many eigenvalues are there in this line
        eigcount = char_count(clearspace(dummy(8:8)//dummy(21:21)//dummy(34:34)//dummy(47:47)//dummy(60:60)))
        if (eigcount.gt.0) then
         
           write(form,"(A1,I1,A6)")'(',eigcount,'F13.7)'

           !write(*,*)dummy(5:70)
           read(dummy(5:70),form)Eigenvalues_tot(eigcount_tot+1:eigcount_tot+eigcount)!Eigenvalues(1:eigcount)
           !Eigenvalues_tot(eigcount_tot+1:eigcount_tot+eigcount) = Eigenvalues
           eigcount_tot = eigcount_tot + eigcount
           !write(*,*)Eigenvalues_tot(1:eigcount_tot)
           !write(*,*)eigenvalues
        endif
     endif
     if (dummy(6:20).eq."EIGENVALUES ARE") then
        !turn the eigenvalue read on
        kidx = kidx + 1
        eigcount_tot = 0
        Eigenvalues_tot(:) = 0d0
        kread = .true.
        !write(*,*)"kidx",kidx
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
   !write(*,*)Eigenvalues_tot(j1)
   Eigenvalues_tot(j1) = Eigenvalues_tot(j1) - Fermienergy
   Eigenvalues_tot(j1) = Eigenvalues_tot(j1)*13.6d0
   !write(*,*)Eigenvalues_tot(j1),Emin,Emax
   if ((Eigenvalues_tot(j1).ge.Emin).and.(Eigenvalues_tot(j1).le.Emax).and.(flagmin.eqv..false.))then
      minidx = j1
      flagmin = .true.

   endif
  if (Eigenvalues_tot(j1).le.Emax) then
      maxidx = j1
   endif
   
 enddo
 if (flagmin.eqv..false.) then
    !no bands in window for this k point
    write(*,"(A2,I4,A9,I3,A18,A43)")"k=",kidx," | #eig.=",eigcount_tot," | min in int.= --",&
                                     " | max in int.= -- | #bands in en.-int.=  0"
 else
    write(*,"(A2,I4,A9,I3,A15,I3,A15,I3,A22,I3)")"k=",kidx," | #eig.=",eigcount_tot," | min in int.=",minidx,&
                                     " | max in int.=",maxidx," | #bands in en.-int.=",maxidx-minidx+1
 endif

END SUBROUTINE
