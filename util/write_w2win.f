PROGRAM write_w2win
 !program prepares input file seed.w2win for wien2wannier 
 ! P.Wissgott 10/01/10
 use util 

 implicit none
 character*50 seedname

  integer iarg,argcount !command line input argument counter
  integer :: nemin,nemax
  integer :: Ljmax  !plain wave expansion cutoff
  integer :: nproj  !#projection = # final Wannier functions
  integer :: scounter,pcounter,dcounter,fcounter                !orbital counters
  integer :: norbitals_s,norbitals_p(3),norbitals_d(5),norbitals_f(7) !#entries per orbital
  integer :: mindex_s,mindex_p(3,2),mindex_d(5,2),mindex_f(7,2)      !m quantum number for orbitals
  integer :: j1,j2,j            !counter
  integer :: atomidx
  character*1 orbital_character
  complex*8 :: slist(1,2)  !entry coefficients
  complex*8 :: plist(3,2)  !entry coefficients
  complex*8 :: dlist(5,2)  !entry coefficients
  complex*8 :: flist(7,2)  !entry coefficients
  character*50  :: sprojtext
  character*50, dimension(1:3) :: pprojtext
  character*50, dimension(1:5) :: dprojtext
  character*50, dimension(1:7) :: fprojtext
  type(structure) lattice
  character*70 dummy,startmessage
  character*8 w2winfileend
  
  !default fileending: non spin-polarized
  w2winfileend = ".w2win"
  startmessage = "++ Preparing a non spin-polarized input file ++"

 !command line argument read-in
 iarg=iargc()
 argcount = 1
 if(iarg.ge.1) then
     do j=1,iarg
        call getarg(j,dummy)
        if (dummy(1:1).eq.'-') then
           if ((dummy(2:3).eq.'up').or.(dummy(2:3).eq.'dn')) then     
              !for spin-polarized calc. the fileendings have additional up/dn
              w2winfileend = ".w2win"//dummy(2:3)
              startmessage = "++ Preparing a spin-polarized input file:"//dummy(2:3)//" ++"
           elseif (dummy(2:2).eq.'h') then
               write(*,*) 'prepares input file case.w2win[up/dn] for wien2wannier'
               write(*,*) 'Usage: write_w2win [-up/-dn] case'
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
     write(*,*) 'Error: At least "seed" commandline argument has to be given.'
     stop
  endif

  write(*,*)startmessage
 
 
 !create the lists for the orbitals
 !s orbital
 norbitals_s = 1
 mindex_s = 0
 slist(1,1) = cmplx(1d0,0d0)
 slist(1,2) = cmplx(0d0,0d0)
 scounter = 1
 sprojtext = "s orbital"

 !p orbitals
 norbitals_p = (/2,2,1/)
 mindex_p(1,:) = (/ -1, 1/)
 mindex_p(2,:) = (/ -1, 1/)
 mindex_p(3,:) = (/  0, 0/)
 plist(1,1) = cmplx(1d0/sqrt(2d0),0d0)!px
 plist(1,2) = cmplx(-1d0/sqrt(2d0),0d0)!px
 pprojtext(1) = "p-x orbital"
 plist(2,1) = cmplx(0d0,1d0/sqrt(2d0))!py
 pprojtext(2) = "p-y orbital"
 plist(2,2) = cmplx(0d0,1d0/sqrt(2d0))!py
 plist(3,1) = cmplx(1d0,0d0)!pz
 plist(3,2) = cmplx(0d0,0d0)!pz
 pprojtext(3) = "p-z orbital"
 pcounter = 1

 !d orbitals
 norbitals_d = (/2,2,2,2,1/)
 mindex_d(1,:) = (/-2, 2/)
 mindex_d(2,:) = (/ -1, 1/)
 mindex_d(3,:) = (/ -1, 1/)
 mindex_d(4,:) = (/-2, 2/)
 mindex_d(5,:) = (/0, 0/)
 dlist(1,1) = cmplx(0d0,1d0/sqrt(2d0))!dxy
 dlist(1,2) = cmplx(0d0,-1d0/sqrt(2d0))!dxy
 dprojtext(1) = "d-xy orbital"
 dlist(2,1) = cmplx(0d0,1d0/sqrt(2d0))!dyz
 dlist(2,2) = cmplx(0d0,1d0/sqrt(2d0))!dyz
 dprojtext(2) = "d-yz orbital"
 dlist(3,1) = cmplx(1d0/sqrt(2d0),0d0)!dxz
 dlist(3,2) = cmplx(-1d0/sqrt(2d0),0d0)!dxz
 dprojtext(3) = "d-xz orbital"
 dlist(4,1) = cmplx(1d0/sqrt(2d0),0d0)!dx2-y2
 dlist(4,2) = cmplx(1d0/sqrt(2d0),0d0)!dx2-y2
 dprojtext(4) = "d-(x2-y2) orbital"
 dlist(5,1) = cmplx(1d0,0d0)!dz2
 dlist(5,2) = cmplx(0d0,0d0)!dz2
 dprojtext(5) = "d-z2 orbital"
 dcounter = 1
 

 !f orbitals
 norbitals_f = (/2,2,2,2,2,2,1/)
 mindex_f(1,:) = (/-3, 3/)
 mindex_f(2,:) = (/-2, 2/)
 mindex_f(3,:) = (/ -1, 1/)
 mindex_f(4,:) = (/ -1, 1/)
 mindex_f(5,:) = (/-2, 2/)
 mindex_f(6,:) = (/-3, 3/)
 mindex_f(7,:) = (/0, 0/)
 flist(1,1) = cmplx(0d0,1d0/sqrt(2d0))!fy(3x2-y2))
 flist(1,2) = cmplx(0d0,1d0/sqrt(2d0))!fy(3x2-y2))
 fprojtext(1) = "f-y(3x2-y2) orbital"
 flist(2,1) = cmplx(1d0/sqrt(2d0),0d0)!fz(x2-3y2))
 flist(2,2) = cmplx(1d0/sqrt(2d0),0d0)!fz(x2-3y2))
 fprojtext(2) = "f-z(x2-y2) orbital"
 flist(3,1) = cmplx(0d0,1d0/sqrt(2d0))!fyz2
 flist(3,2) = cmplx(0d0,1d0/sqrt(2d0))!fyz2
 fprojtext(3) = "f-yz2 orbital"
 flist(4,1) = cmplx(1d0/sqrt(2d0),0d0)!fxz2
 flist(4,2) = cmplx(-1d0/sqrt(2d0),0d0)!fxz2
 fprojtext(4) = "f-xz2 orbital"
 flist(5,1) = cmplx(0d0,1d0/sqrt(2d0))!fxyz
 flist(5,2) = cmplx(0d0,-1d0/sqrt(2d0))!fxyz
 fprojtext(5) = "f-xyz orbital"
 flist(6,1) = cmplx(1d0/sqrt(2d0),0d0)!fx(x2-3y2))
 flist(6,2) = cmplx(-1d0/sqrt(2d0),0d0)!fx(x2-3y2))
 fprojtext(6) = "f-x(x2-3y2)) orbital"
 flist(7,1) = cmplx(1d0,0d0)!fz3
 flist(7,2) = cmplx(0d0,0d0)!fz3
 fprojtext(7) = "f-z3 orbital"
 fcounter = 1  

 !read-in from seed.struct/seed.outputkgen
 open(unit=13,file=clearspace(seedname)//'.struct',status='old')
 open(unit=14,file=clearspace(seedname)//'.outputkgen',status='old')
 call countatoms(13,lattice)
 call readin_lattice(13,14,lattice)
 close(13)
 close(14)

 open(unit=11,file=clearspace(seedname)//clearspace(w2winfileend),status='unknown')

 Ljmax = 3
 write(11,"(A4)")"BOTH"
 
 print*,"Enter minimal and maximal Bloch band,[n1 n2]:"
 read (*,*) nemin,nemax
 if (nemax.lt.nemin) then
   write(*,*)"Error: minimal band number has to be smaller/larger than maximal"
   stop
 endif
 write(11,"(2I4,A60)")nemin,nemax,      "   # min band Nmin, max band Nmax                         "

 print*,"Enter number of Wannier functions,[n1]:"
 read (*,*) nproj
 if (nproj.gt.nemax-nemin+1) then
   write(*,*)"Error: number of Wannier function must not be larger than number of Bloch bands"
   stop
 endif 
 write(11,"(2I3,A60)")Ljmax,nproj,"     # LJMAX max in exp(ibr) expansion, #Wannier functions"

 write(*,*)"Setting initial values for Wannier functions:"
 do j1=1,nproj
    write(*,"(A8,I2)")"Orbital ",j1
    print*,"Enter center atom and character[n1 s/p/d/f]:(e.g.: 2 d)"
    read (*,*) atomidx,orbital_character
    if (atomidx.gt.lattice%atomcount) then
       write(*,*)"Error: atom index out of bounds"
       stop
    endif
    select case (orbital_character)
       case ("s")
          
          write(11,"(I1,A60)")norbitals_s,"     #"//sprojtext
          do j2=1,norbitals_s
             write(11,"(2I2,I3,2F13.8,A47)")atomidx,0,mindex_s,real(slist(scounter,j2)),aimag(slist(scounter,j2)),&
                                              " # index of atom, L, M, coefficient (complex)"
          enddo
          write(*,*)">> added a "//sprojtext
     
       case ("p")
          write(11,"(I1,A60)")norbitals_p(pcounter),"     #"//pprojtext(pcounter)
          do j2=1,norbitals_p(pcounter)
             write(11,"(2I2,I3,2F13.8,A47)")atomidx,1,mindex_p(pcounter,j2),real(plist(pcounter,j2)),aimag(plist(pcounter,j2)),&
                                              " # index of atom, L, M, coefficient (complex)"
          enddo
          write(*,*)">> added a "//trim(pprojtext(pcounter))//"(local coordinate system)"
       
          if (pcounter.eq.3) then 
             pcounter = 1
          else
             pcounter = pcounter + 1
          endif

       case ("d")
          write(11,"(I1,A60)")norbitals_d(dcounter),"     #"//dprojtext(dcounter)
          do j2=1,norbitals_d(dcounter)
             write(11,"(2I2,I3,2F13.8,A47)")atomidx,2,mindex_d(dcounter,j2),real(dlist(dcounter,j2)),aimag(dlist(dcounter,j2)),&
                                              " # index of atom, L, M, coefficient (complex)"
          enddo
          write(*,*)">> added a "//trim(dprojtext(dcounter))//"(local coordinate system)"

          if (dcounter.eq.5) then 
             dcounter = 1
          else
             dcounter = dcounter + 1
          endif

        case ("f")
          write(11,"(I1,A60)")norbitals_f(fcounter),"     #"//fprojtext(fcounter)
          do j2=1,norbitals_f(fcounter)
             write(11,"(2I2,I3,2F13.8,A47)")atomidx,3,mindex_f(fcounter,j2),real(flist(fcounter,j2)),aimag(flist(fcounter,j2)),&
                                              " # index of atom, L, M, coefficient (complex)"
          enddo
          write(*,*)">> added a "//trim(fprojtext(fcounter))//"(local coordinate system)"
       
          if (fcounter.eq.7) then 
             fcounter = 1
          else
             fcounter = fcounter + 1
          endif

   case default
      write(*,*)"Error: Unknown orbital character"
      stop
   end select

 enddo
 write(11,*)
 close(11)

END PROGRAM



