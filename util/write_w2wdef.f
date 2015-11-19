PROGRAM write_w2wdef
 !program prepares input w2w.def for wien2wannier 
 ! P.Wissgott 10/01/10

 use util 

 implicit none
 character*50 seedname

 integer iarg !command line input argument counter
 integer :: j,argcount !counter
 character*70 argdummy,startmessage
 character*8 w2winfileend,scfwffileend
 character*9 w2woutfileend
 character*11 energyfileend,vectorfileend
 character*6 amnfileend,mmnfileend,eigfileend,vspfileend,deffileend

 !default fileending: non spin-polarized
  w2winfileend = ".w2win"
  vectorfileend = ".vector"
  scfwffileend  = ".scfwf"
  w2woutfileend = ".w2wout"
  energyfileend = ".energy"
  amnfileend = ".amn"  
  mmnfileend = ".mmn"
  eigfileend = ".eig"
  vspfileend = ".vsp"
  deffileend = ".def"
  startmessage = "++ Preparing a non spin-polarized def file ++"
 
 !command line argument read-in
 iarg=iargc()
 argcount = 1
 if(iarg.ge.1) then
     do j=1,iarg
        call getarg(j,argdummy)
        if (argdummy(1:1).eq.'-') then
           if ((argdummy(2:3).eq.'up').or.(argdummy(2:3).eq.'dn')) then     
              !for spin-polarized calc. the fileendings have additional up/dn
              w2winfileend = ".w2win"//argdummy(2:3)
              vectorfileend = ".vector"//argdummy(2:3)
              scfwffileend  = ".scfwf"//argdummy(2:3)
              w2woutfileend = ".w2wout"//argdummy(2:3)
              energyfileend = ".energy"//argdummy(2:3)
              amnfileend = ".amn"//argdummy(2:3)
              mmnfileend = ".mmn"//argdummy(2:3)
              eigfileend = ".eig"//argdummy(2:3)
              vspfileend = ".vsp"//argdummy(2:3)
              startmessage = "++ Preparing a spin-polarized def file:"//argdummy(2:3)//" ++"
           elseif ((argdummy(2:5).eq.'soup').or.(argdummy(2:5).eq.'sodn')) then     
              !for SO calc. 
              w2winfileend = ".w2win"//argdummy(4:5)
              vectorfileend = ".vector"//argdummy(2:5)
              scfwffileend  = ".scfwf"//argdummy(4:5)
              w2woutfileend = ".w2wout"//argdummy(4:5)
              energyfileend = ".energy"//argdummy(2:5)
              amnfileend = ".amn"//argdummy(4:5)
              mmnfileend = ".mmn"//argdummy(4:5)
              eigfileend = ".eig"//argdummy(4:5)
              vspfileend = ".vsp"//argdummy(4:5)
              startmessage = "++ Preparing a spin-polarized def file for SO:"//argdummy(4:5)//" ++"
         elseif (argdummy(2:2).eq.'h') then
               write(*,*) 'prepares input file w2w.def for wien2wannier'
               write(*,*) 'Usage: write_w2wdef [-up/-dn/-soup/-sodn] case'
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
 
 open(unit=1,file='w2w'//clearspace(deffileend),status='unknown')

 write(1,*)" 5,'",clearspace(seedname),clearspace(w2winfileend),"',   'old',    'formatted',0"
 write(1,*)" 6,'",clearspace(seedname),clearspace(w2woutfileend),"','unknown','formatted',0"
 write(1,*)" 7,'",clearspace(seedname),clearspace(amnfileend),"','unknown','formatted',0"
 write(1,*)" 8,'",clearspace(seedname),clearspace(mmnfileend),"','unknown','formatted',0"
 write(1,*)" 9,'",clearspace(seedname),clearspace(vectorfileend),"','unknown','unformatted',9000"
 write(1,*)"10,'",clearspace(seedname),".nnkp','old','formatted',0"
 write(1,*)"12,'",clearspace(seedname),clearspace(eigfileend),"','unknown','formatted',0"
 write(1,*)"18,'",clearspace(seedname),clearspace(vspfileend),"',     'old','formatted',0"
 write(1,*)"20,'",clearspace(seedname),".struct',   'old',    'formatted',0"
 write(1,*)"21,'",clearspace(seedname),clearspace(scfwffileend),"',   'unknown','formatted',0"
 write(1,*)"50,'",clearspace(seedname),clearspace(energyfileend),"', 'unknown','formatted',9000"

 close(1)

END PROGRAM



