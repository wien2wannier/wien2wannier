PROGRAM write_w2wdef
 !program prepares input w2w.def for wien2wannier 
 ! P.Wissgott 10/01/10
 
 use util 

 implicit none
 character*50 seedname
 integer iarg !command line input argument counter
 integer argcount,j
 character*70 argdummy,startmessage
 character*8 w2winfileend,psinkfileend
 character*6 vspfileend,abcfileend,chkfileend
 character*7 gridfileend
 character*10 wplotinfileend
 character*11 wplotoutfileend,vectorfileend
 character*9 psiargfileend


 w2winfileend = ".w2win"
 vspfileend = ".vsp"
 gridfileend  = ".grid"
 wplotinfileend = ".wplotin"
 wplotoutfileend = ".wplotout"
 abcfileend = ".abc"
 psinkfileend = ".psink"
 chkfileend = ".chk"
 psiargfileend = ".psiarg"
 vectorfileend = ".vector"
 startmessage = "++ Preparing a standard def file ++" 

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
              vspfileend = ".vsp"//argdummy(2:3)
              gridfileend  = ".grid"//argdummy(2:3)
              wplotinfileend = ".wplotin"//argdummy(2:3)
              wplotoutfileend = ".wplotout"//argdummy(2:3)
              abcfileend = ".abc"//argdummy(2:3)
              psinkfileend = ".psink"//argdummy(2:3)
              chkfileend = ".chk"//argdummy(2:3)
              psiargfileend = ".psiarg"//argdummy(2:3)
              vectorfileend = ".vector"//argdummy(2:3)
              startmessage = "++ Preparing a spin-polarized def file:"//argdummy(2:3)//" ++"
          elseif (argdummy(2:5).eq.'soup') then     
              !for SO calc. 
              w2winfileend = ".w2winup"
              vspfileend = ".vspup"
              gridfileend  = ".gridup"
              wplotinfileend = ".wplotinup"
              wplotoutfileend = ".wplotoutup"
              abcfileend = ".abcup"
              psinkfileend = ".psinkup"
              chkfileend = ".chkso"
              psiargfileend = ".psiargup"
              vectorfileend = ".vectorsoup"
              startmessage = "++ Preparing a spin-polarized SO def file: soup ++"
          elseif (argdummy(2:5).eq.'sodn') then     
              !for SO calc. 
              w2winfileend = ".w2windn"
              vspfileend = ".vspdn"
              gridfileend  = ".griddn"
              wplotinfileend = ".wplotindn"
              wplotoutfileend = ".wplotoutdn"
              abcfileend = ".abcdn"
              psinkfileend = ".psinkdn"
              chkfileend = ".chkso"
              psiargfileend = ".psiargdn"
              vectorfileend = ".vectorsodn"
              startmessage = "++ Preparing a spin-polarized SO def file: sodn ++"
           elseif (argdummy(2:2).eq.'h') then
               write(*,*) 'prepares input file wplot.def for wplot'
               write(*,*) 'Usage: write_wplotdef [-up/-dn/-soup/-sodn] case'
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

 open(unit=1,file='wplot.def',status='unknown')

 write(1,*)" 2,'",clearspace(seedname),".tmp',   'unknown',    'formatted',0"
 write(1,*)" 5,'",clearspace(seedname),clearspace(wplotinfileend),"',   'old','formatted',0"
 write(1,*)" 6,'",clearspace(seedname),clearspace(wplotoutfileend),"','unknown','formatted',0"
 write(1,*)" 7,'",clearspace(seedname),clearspace(gridfileend),"','unknown','formatted',0"
 write(1,*)" 8,'",clearspace(seedname),".struct',   'old',    'formatted',0" 
 write(1,*)"10,'",clearspace(seedname),clearspace(vectorfileend),"','old','unformatted',0"
 write(1,*)"12,'",clearspace(seedname),clearspace(abcfileend),"','unknown','unformatted',0"
 write(1,*)"18,'",clearspace(seedname),clearspace(vspfileend),"','old','formatted',0"
 write(1,*)"21,'",clearspace(seedname),clearspace(psinkfileend),"',   'unknown','formatted',0"
 write(1,*)"22,'",clearspace(seedname),clearspace(psiargfileend),"', 'unknown','formatted',0"
 write(1,*)"31,'",clearspace(seedname),clearspace(w2winfileend),"',   'unknown','formatted',0"
 write(1,*)"32,'",clearspace(seedname),clearspace(chkfileend),"', 'old','unformatted',0"
 
 close(1)

END PROGRAM



