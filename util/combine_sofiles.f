PROGRAM write_w2win
 !program prepares input file seed.w2win for wien2wannier 
 
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

 !p orbitals
 norbitals_p = (/2,2,1/)
 mindex_p(1,:) = (/ -1, 1/)
 mindex_p(2,:) = (/ -1, 1/)
 mindex_p(3,:) = (/  0, 0/)
 plist(1,1) = cmplx(0d0,1d0/sqrt(2d0))
 plist(1,2) = cmplx(0d0,-1d0/sqrt(2d0))
 plist(2,1) = cmplx(1d0/sqrt(2d0),0d0)
 plist(2,2) = cmplx(1d0/sqrt(2d0),0d0)
 plist(3,1) = cmplx(1d0,0d0)
 plist(3,2) = cmplx(0d0,0d0)
 pcounter = 1

 !d orbitals
 norbitals_d = (/2,2,2,2,1/)
 mindex_d(1,:) = (/-2, 2/)
 mindex_d(2,:) = (/ -1, 1/)
 mindex_d(3,:) = (/ -1, 1/)
 mindex_d(4,:) = (/-2, 2/)
 mindex_d(5,:) = (/0, 0/)
 dlist(1,1) = cmplx(0d0,1d0/sqrt(2d0))
 dlist(1,2) = cmplx(0d0,-1d0/sqrt(2d0))
 dlist(2,1) = cmplx(1d0/sqrt(2d0),0d0)
 dlist(2,2) = cmplx(1d0/sqrt(2d0),0d0)
 dlist(3,1) = cmplx(0d0,1d0/sqrt(2d0))
 dlist(3,2) = cmplx(0d0,-1d0/sqrt(2d0))
 dlist(4,1) = cmplx(1d0/sqrt(2d0),0d0)
 dlist(4,2) = cmplx(1d0/sqrt(2d0),0d0)
 dlist(5,1) = cmplx(1d0,0d0)
 dlist(5,2) = cmplx(0d0,0d0)
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
 flist(1,1) = cmplx(0d0,1d0/sqrt(2d0))
 flist(1,2) = cmplx(0d0,-1d0/sqrt(2d0))
 flist(2,1) = cmplx(1d0/sqrt(2d0),0d0)
 flist(2,2) = cmplx(1d0/sqrt(2d0),0d0)
 flist(3,1) = cmplx(0d0,1d0/sqrt(2d0))
 flist(3,2) = cmplx(0d0,-1d0/sqrt(2d0))
 flist(4,1) = cmplx(0d0,1d0/sqrt(2d0))
 flist(4,2) = cmplx(0d0,-1d0/sqrt(2d0))
 flist(5,1) = cmplx(1d0/sqrt(2d0),0d0)
 flist(5,2) = cmplx(1d0/sqrt(2d0),0d0)
 flist(6,1) = cmplx(0d0,1d0/sqrt(2d0))
 flist(6,2) = cmplx(0d0,-1d0/sqrt(2d0))
 flist(7,1) = cmplx(1d0,0d0)
 flist(7,2) = cmplx(0d0,0d0)
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
          
          write(11,"(I1,A60)")norbitals_s,"     #number of entries per orbital (list of projected orbitals)"
          do j2=1,norbitals_s
             write(11,"(2I2,I3,2F13.8,A47)")atomidx,0,mindex_s,real(slist(scounter,j2)),aimag(slist(scounter,j2)),&
                                              " # index of atom, L, M, coefficient (complex)"
          enddo
     
       case ("p")
          write(11,"(I1,A60)")norbitals_p(pcounter),"     #number of entries per orbital (list of projected orbitals)"
          do j2=1,norbitals_p(pcounter)
             write(11,"(2I2,I3,2F13.8,A47)")atomidx,1,mindex_p(pcounter,j2),real(plist(pcounter,j2)),aimag(plist(pcounter,j2)),&
                                              " # index of atom, L, M, coefficient (complex)"
          enddo
       
          if (pcounter.eq.3) then 
             pcounter = 1
          else
             pcounter = pcounter + 1
          endif

       case ("d")
          write(11,"(I1,A60)")norbitals_d(dcounter),"     #number of entries per orbital (list of projected orbitals)"
          do j2=1,norbitals_d(dcounter)
             write(11,"(2I2,I3,2F13.8,A47)")atomidx,2,mindex_d(dcounter,j2),real(dlist(dcounter,j2)),aimag(dlist(dcounter,j2)),&
                                              " # index of atom, L, M, coefficient (complex)"
          enddo
       
          if (dcounter.eq.5) then 
             dcounter = 1
          else
             dcounter = dcounter + 1
          endif

        case ("f")
          write(11,"(I1,A60)")norbitals_f(fcounter),"     #number of entries per orbital (list of projected orbitals)"
          do j2=1,norbitals_f(fcounter)
             write(11,"(2I2,I3,2F13.8,A47)")atomidx,3,mindex_f(fcounter,j2),real(flist(fcounter,j2)),aimag(flist(fcounter,j2)),&
                                              " # index of atom, L, M, coefficient (complex)"
          enddo
       
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
 close(11)

END PROGRAM



