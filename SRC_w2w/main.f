!!! wien2wannier/SRC_w2w/main.f
!!!
!!!    Main program ‘w2w’
!!!
!!! Copyright 2010-2012 Jan Kuneš, Philipp Wissgott
!!!           2013-2015 Elias Assmann

!!/--- Files expected in ‘def’ ---
!!  5 inwf	'old'	  'formatted'
!!  6 outputwf	'unknown' 'formatted'
!!  7 amn	'unknown' 'formatted'
!!  8 mmn	'unknown' 'formatted'
!!  9 vector	'unknown' 'unformatted'
!! 10 nnkp	'old'	  'formatted'
!! 12 eig	'unknown' 'formatted'
!! 18 vsp	'old'	  'formatted'
!! 20 struct	'old'	  'formatted'
!! 50 energy	'old'     'formatted'
!! 51 fermi	'old'     'formatted'
!!\---

program wf
  use param
  use const,    only: R8, BUFSZ
  use struct,   only: aa,bb,cc, irel, alpha, Nat, lattic, title, init_struct
  use xa,       only: init_xa
  use xa3,      only: init_xa3
  use bessel,   only: init_bessel
  use Amn_Mmn,  only: c, lmax2, Nmat, Nrad, Nrf, init_Mmn
  use ams,      only: init_ams
  use pairs,    only: kp, kpb, bqx,bqy,bqz, bqx1,bqy1,bqz1, init_pairs
  use util,     only: paropen
  use wien2k,   only: errflg, errclr, gtfnam
  use gener,    only: br1, br2

  implicit none

  character(len=    3)  :: mode
  character(len=   11)  :: status,form
  character(len=   67)  :: errmsg
  character(len=BUFSZ)  :: deffn, errfn, aline
  character(len=BUFSZ)  :: fname, vecfn, enefn

  logical :: Mmn, Amn
  integer :: centeratom(300)
  integer :: iloop, i, ia, ii, ind, info, ios, iproc, irecl, iunit, j, l, m
  integer :: ljmax, nproj, maxx,maxy,maxz, maxg, n, n_pair, nen
  integer :: nemin, nemax, Nb, Nk, nntot, bx,by,bz, kkk

  real(r8) :: t1, t2, t3, x1, x2
  real(r8) :: efermi

!-----------------------------------------------------------------------

  call init_ams
  call gtfnam(deffn,errfn,iproc)
  call errflg(errfn,'Error in W2W')
  open(unit_def, FILE=deffn, STATUS='old', ERR=910)
  def: do
     read(unit_def, *, END=20, ERR=960) iunit, fname, status, form, irecl

     select case (iunit)
     case (unit_vector)
        vecfn=fname
     case (50)
        enefn=fname
     case default
        open(iunit, FILE=fname, STATUS=status, FORM=form, ERR=920)
     end select
  end do def
20 close(unit_def)

  write(unit_out, '("W2W ", A /)'), wien2wannier_version

!!!.....READ STRUCT
  CALL init_struct
!!!....Find nmat and Nk
!!! Nk could be gotten easier from ‘klist’ -- can we get nmat
!!! somewhere else?  (‘vector’?)
  Nk=0

  enefile: do iloop=1,max(iproc,1)
     call paropen(unit_ene, ENEFN, iproc, iloop, STATUS='old')

     header: DO I=1,NAT
        READ(unit_ene,*)
        READ(unit_ene,*)
     END DO header
     kpts: DO
        READ(unit_ene,'(67X, 2i6)', IOSTAT=ios) N, NEn
        IF (ios /= 0) EXIT kpts
        Nk=Nk+1
        nmat=MAX(n,nmat)
        DO ii=1,nen
           READ(unit_ene,*)
        ENDDO
     END DO kpts

     CLOSE(unit_ene)
  end do enefile

!!!.....READ INPUT AND POTE
  nnkpt: do
     read(unit_nnkp,'(a80)',end=114) aline
     if(index(aline,'begin kpoints').ne.0) then
        read(unit_nnkp,*,err=913) n
        if (n/=Nk) goto 941
        cycle nnkpt
     endif

     if(index(aline,'begin nnkpts').ne.0) then
        read(unit_nnkp,*,err=914) nntot
        n_pair=Nk*nntot
        call init_pairs(n_pair)
        do i=1,N_pair
           read(unit_nnkp,*) kp(i),kpb(i),bqx1(i),bqy1(i),bqz1(i)
        enddo
        exit nnkpt
     endif
  end do nnkpt
114 close(unit_nnkp)

  write(unit_out,*)'NUM_KPTS=',Nk
  write(unit_out,*)'NNTOT=',NNTOT
  write(unit_out,*)'N_pair=',N_pair

  MMN=.true.
  AMN=.true.
  READ(unit_in,*)MODE
  IF (MODE.eq.'MMN') AMN=.false.
  IF (MODE.eq.'AMN') MMN=.false.
  write(unit_out,*)'MODE='
  if (MMN) write(unit_out,*)'MMN'
  if (AMN) write(unit_out,*)'AMN'
  READ(unit_in,*)NEMIN,NEMAX
  if (nemax.lt.nemin) stop 'nemin > nemax'
  READ(unit_in,*)LJMAX,NPROJ
  write(unit_out,*)'nemin,nemax:',nemin,nemax
  Nb=nemax-nemin+1
  CALL init_xa3(Nb,nmat,Nk+1)
  CALL init_xa(LMAX2,NMAT,NRAD,Nb)
  CALL init_mmn(Nb,n_pair,nproj)

  read(unit_fermi, *) efermi

  kkk=0
  maxx=0; maxy=0; maxz=0
  vectorfiles: do iloop=1, max(iproc, 1)
     call paropen(unit_vector, vecfn, iproc, iloop, &
          STATUS='old', FORM='unformatted')

     CALL read_vec(nemin,nemax,kkk,maxx,maxy,maxz,efermi)

     CLOSE(unit_vector)
  end do vectorfiles

  if(kkk /= Nk) goto 942

  C(:,:,:) = dcmplx(0d0,0d0)

  read_proj: IF (AMN) THEN
     DO I=1,NPROJ
        READ(unit_in,*,err=925,end=926)N
        DO J=1,N
           READ(unit_in,*,err=925,end=926)IA,L,M,X1,X2
           ind=L*(L+1)+M+1
           C(I,ind,IA)=X1 + (0,1)*X2
           CENTERATOM(I)=IA
        ENDDO
     ENDDO

     write(unit_out,*)'initial orbital projections'
     DO I=1,nproj
        write(unit_out,*)'orbital #',I,'centered at atom',CENTERATOM(I)
        ind=0
        DO L=0,3
           DO M=-L,L
              ind=ind+1
              write(unit_out,*)L,M,C(I,ind,CENTERATOM(I))
           ENDDO
        ENDDO
     ENDDO
  ENDIF read_proj

  call init_bessel(LMAX2,LJMAX,NRAD,NRF)
  call gaunt2
  WRITE(unit_out,800)
  WRITE(unit_out,805)  TITLE
  if (amn) then
     write(unit_amn,806) TITLE
     write(unit_amn,807) NEMAX-NEMIN+1,Nk,NPROJ
  endif
  if (mmn) then
     write(unit_mmn,806) TITLE
     write(unit_mmn,807) NEMAX-NEMIN+1,Nk,NNTOT
  endif

  WRITE(unit_out,810)  LATTIC
  WRITE(unit_out,820)  AA,BB,CC
  WRITE(unit_out,840)  NAT
  WRITE(unit_out,850)  IREL

  CALL LATGEN
  !     rotate boundary vectors
  do i=1,N_pair
     bx=bqx1(i); by=bqy1(i); bz=bqz1(i)
     bqx(i) = int(bx*br2(1,1) + by*br2(1,2) + bz*br2(1,3))
     bqy(i) = int(bx*br2(2,1) + by*br2(2,2) + bz*br2(2,3))
     bqz(i) = int(bx*br2(3,1) + by*br2(3,2) + bz*br2(3,3))
  enddo
  write(unit_out,*)' alpha test',(alpha(i),i=1,3)
  !.....CALCULATE CHARGE DENSITY CLM(R) IN SPHERES,  PARTIAL CHARGES

  ! l2mmn, l2amn need vector file to be open
  call paropen(unit_vector, vecfn, iproc, 1, STATUS='old', FORM='unformatted')

  if (MMN) then
     call cputim(t1)
     call l2mmn(Nb,Nk,NNTOT,LJMAX)
     call cputim(t2)
     !JXZ: MAXG is not pre-defined
     MAXG = 0
     write(unit_out,*)'MXG=',MAXG
     call planew(Nb,Nk,NNTOT,maxx+1,maxy+1,maxz+1)
     call cputim(t3)
     write(unit_out,*)'CPU l2mmn:',t2-t1
     write(unit_out,*)'CPU planew:',t3-t2
  endif

  if (AMN) then
     call cputim(t1)
     call l2amn(Nb,NPROJ,Nk)
     call cputim(t2)
     write(unit_out,*)'CPU l2amn:',t2-t1
  endif

  CALL ERRCLR(ERRFN)
  STOP 'W2W END'

!!!        error handling
!!!
910 INFO = 1

!!!        def file couldn't be opened
!!!
  WRITE (ERRMSG,9000) DEFFN
  CALL OUTERR('w2w',ERRMSG)
  GOTO 999
920 INFO = 2

!!!        file FNAME couldn't be opened
!!!
  WRITE (ERRMSG,9010) IUNIT
  CALL OUTERR('w2w',ERRMSG)
  WRITE (ERRMSG,9020) FNAME
  CALL OUTERR('w2w',ERRMSG)
  WRITE (ERRMSG,9030) STATUS, FORM
  CALL OUTERR('w2w',ERRMSG)
  GOTO 999
925 WRITE(*,*)'error reading projection',I
  GOTO 999
926 WRITE(*,*)'too few projections'
960 info = 7
913 call outerr('w2w', 'error reading begin kpoints')
  goto 999
914 call outerr('w2w', 'error reading begin nnkpts')
  goto 999
941 call outerr('w2w', 'inconsistent numbers of k-points between nnkp and energy files')
  goto 999
942 call outerr('w2w', 'inconsistent numbers of k-points between vector and energy files')
  goto 999

!!!        Error reading file 'lapw2.def'
!!!
  WRITE (ERRMSG,9040) FNAME
  CALL OUTERR('w2w',ERRMSG)
  GOTO 999
999 STOP 'w2w - Error'
  !
  !
800 FORMAT(////,30X,50(1H-),/,33X,'S T R U C T U R A L   ', &
         & 'I N F O R M A T I O N',/,30X,50(1H-),//)
805 FORMAT(3X,'SUBSTANCE',20X,'= ',A80,/)
806 FORMAT(A20)
807 FORMAT(3I12)
810 FORMAT(3X,'LATTICE',22X,'= ',A4)
820 FORMAT(3X,'LATTICE CONSTANTS ARE',8X,'= ',3F12.7)
840 FORMAT(3X,'NUMBER OF ATOMS IN UNITCELL  = ',I3)
850 FORMAT(3X,'MODE OF CALCULATION IS',7X,'= ',A4)
9000 FORMAT('can''t open definition file ',A40)
9010 FORMAT('can''t open unit: ',I2)
9020 FORMAT('       filename: ',A50)
9030 FORMAT('         status: ',A,'  form: ',A)
9040 FORMAT('Error reading file: ',A47)

END program wf


!!/---
!! Local Variables:
!! mode: f90
!! End:
!!\---
!!
!! Time-stamp: <2016-04-08 17:27:49 assman@faepop36.tu-graz.ac.at>
