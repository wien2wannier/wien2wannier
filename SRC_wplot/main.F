!!! wien2wannier/SRC_wplot/main.F

subroutine main()
  use struct,    only: pos, Znuc, RMT, mult, iatnr
  use radgrd,    only: rm, rnot, dx, jri
  use latt,      only: Vuc, br1, br2, br3, br4
  use sym2,      only: tau, imat, iord
  use lolog,     only: lapw, ilo, Nlo
  use loabc,     only: Alo
  use atspdt,    only: p, dp
  use radfu,     only: rrad
  use bessfu,    only: fj, dfj, rad, irad
  use work,      only: aug
  use grid,      only: rgrid, ireg, ilat, iri, npg
  use const,     only: DPk, R8, C16, BUFSZ, PI
  use util,      only: paropen, uppercase
  use clio,      only: croak
  use Wannier90, only: chk_t, chk_read
  use structmod, only: struct_t, struct_read
  use param,     only: Nrad, Nsym, NRF, NLOAt, LOmax, Lmax7, idx_wann,      &
       unit_in, unit_out, unit_struct, unit_vector, unit_chk, unit_rot,     &
       unit_inwf, unit_psink, unit_psiarg, outfn, vecfn, psinkfn, psiargfn, &
       iproc

  IMPLICIT none

  character(BUFSZ) :: buf, rotfilename
  character( 3)    :: switch, iunit
  character( 4)    :: lattic, GridMode, RelMode
  character( 5)    :: whpsi
  character(10)    :: aname, bname

  LOGICAL    LARGE
  complex(C16) ALMSUM

  type(chk_t)    :: chk

#ifndef _COMPLEX_
  REAL(R8),     allocatable :: CF(:)
#else
  COMPLEX(C16), allocatable :: CF(:)
#endif

  COMPLEX(C16), allocatable ::  COEF(:) !changed by pwissgott

  real(r8),allocatable :: BK(:,:)

  COMPLEX(C16),allocatable :: PSI(:),tot(:)
  COMPLEX(C16) ALM((LMAX7+1)*(LMAX7+1),NRF), Y((LMAX7+1)*(LMAX7+1))
  COMPLEX(C16) PHS
  real(r8),allocatable ::  ROTLOC(:,:,:),atms(:,:)
  integer,allocatable :: IOP(:),nreg(:)
  integer,allocatable ::  KX(:),KY(:),KZ(:)
  complex(c16), allocatable :: uu(:), rotmatrix(:,:)
  real(r8),allocatable :: tmprot(:)
  LOGICAL :: REL,ORTHO,PRIM,ORTH, DEPHAS,rotate=.false.

  integer :: i, ie, ieabs, ierr, iG, ik, il, il1, imax, index, index1, irf
  integer :: iskpt, iseig, iloop
  integer :: j, j2, j3, jatom, jneq, k, kkk, kkkabs, l, latom, lm, m, mu
  integer :: Nat, Nb, Ndat, NE, NEmax, NEmin, Nmatk, Nmat, Nmt, Npw
  integer :: NP(3), Nskip

  real(R8) :: A(3), ROT0(3,3), R(3), FK(3), phsav, prefac
  real(R8) :: alpha, beta, gamma, arg, E, weight

  complex(kind=DPk) :: local_contribution

! ############################
! # PART I : read input file #
! ############################

!!! We need to pre-read up to SWITCH in order to read ISEIG and open
!!! the output file before we call grdgen() and latgen()
  read(unit_in, '(A3)') switch

  select case (uppercase(switch))
  case ('3D')
     nskip = 5
  case ('2D')
     nskip = 4
  case ('1D')
     nskip = 3
  case default
     nskip = 1
  end select

  do i=1,nskip
     read(unit_in, *)
  end do

!     << read in wave function options >>
  read (unit_in,'(A)') SWITCH
  select case (uppercase(switch))
  case ('DEP')
     switch = 'DEP'
     dephas = .true.
  case ('NO')
     switch = 'NO'
     dephas = .false.
  case default
     call croak("unknown postprocessing option "//switch)
  end select

  read (unit_in, *) IUNIT,WHPSI
  select case (uppercase(iunit))
  case ('', 'AU', 'ATU')
     iunit = 'AU'
     prefac = 1
  case ('ANG')
     iunit = 'ANG'
     prefac = 1 / sqrt( 0.529177_DPk**3 )
  case default
     call croak("unknown units option "//iunit)
  end select

  IF(WHPSI.EQ.'     ')WHPSI='LARGE'
  IF(WHPSI.NE.'SMALL'.AND.WHPSI.NE.'LARGE') &
       call croak('ERROR: UNKNOWN REL. COMPONENT OPTION')
  LARGE = WHPSI .EQ. 'LARGE'

!     << read in wave function selection (and file handling) >>
  READ (unit_in,*) ISKPT,ISEIG
  read(unit_inwf,*)
  read(unit_inwf,*)nemin,nemax
  nb=nemax-nemin+1
  allocate(uu(nb))
  if (idx_wann.gt.0) then
     iseig = idx_wann
  endif

!!! Now we know the WF index ISEIG and can open files accordingly
  call wfopen(unit_psink,  psinkfn,  iseig, STATUS='unknown')
  call wfopen(unit_psiarg, psiargfn, iseig, STATUS='unknown')
  call wfopen(unit_out,    outfn,    iseig, STATUS='unknown')

  write(unit_out,"( ' selected ',A7,' :',I5)") 'Wannier function   ', iseig

!!! Read chk file for num_wann, u_matrix, and friends
  call chk_read(unit_chk, chk)

  write(unit_out,'("Read Wannier90 ‘chk’ file ", A)') trim(chk%header)

  if (chk%num_bands /= nb) &
       call croak("Error in read-in of chk file: &
       &inconsistency of number of Bloch bands")

  if (iseig > chk%num_wann) &
       call croak("Error in read-in of chk file: &
          &number of Wannier functions is smaller than given index")

  close(unit_chk)

!!! Apply WF rotation
  if (rotate) then
     read(unit_in, "(A)") rotfilename
     write(*,*) "rotate with matrix from ", trim(rotfilename)
     open(unit_rot, file=trim(rotfilename), status='old')
     allocate(rotmatrix(chk%num_wann, chk%num_wann), tmprot(chk%num_wann*2))

     do j=1,chk%num_wann
        read(unit_rot,*) tmprot

        do i=1,chk%num_wann
           rotmatrix(j,i) = cmplx(tmprot(2*i-1), tmprot(2*i), DPk)
        enddo
        write(*,"(100F13.5)") rotmatrix(j,:)
     enddo
     close(unit_rot)
  endif

  rewind(unit_in)               ! for grdgen()
  close(unit_inwf)

! ###################################
! # PART II : structural processing #
! ###################################
!
!     << set up the lattice >>
  READ(unit_struct, '(A)') buf
  READ(unit_struct, '(A4,23X,I3/13X,A4)') LATTIC,NAT,RelMode

  allocate ( POS(3,NAT*48),ZNUC(NAT),RMT(NAT),MULT(NAT),IATNR(NAT*48))
  allocate ( RM(NRAD,NAT),RNOT(NAT),DX(NAT),JRI(NAT) )
  allocate ( ILO(0:LOMAX, NAT),LAPW(0:LMAX7,NAT) )
  allocate ( ALO(0:LOMAX,NLOAT,NRF,NAT))
  allocate ( P(0:LMAX7,NRF,NAT),DP(0:LMAX7,NRF,NAT))
  allocate ( RRAD(NRAD,0:LMAX7,NRF,NAT))
  allocate ( RAD(NAT),IRAD(NAT))
  allocate ( ROTLOC(3,3,NAT*48),atms(3,NAT))
  allocate ( IOP(NAT*48),nreg(0:NAT*48))

  REL = RelMode .EQ. 'RELA'
  READ(unit_struct, '(6F10.5)') A,ALPHA,BETA,GAMMA
  IF(ALPHA.EQ.0.0D0) ALPHA=90.0D0
  IF(BETA .EQ.0.0D0) BETA =90.0D0
  IF(GAMMA.EQ.0.0D0) GAMMA=90.0D0

  WRITE(unit_out, '(/A/)') trim(buf)
  WRITE(unit_out,"(1X,A4,' LATTIC WITH ',I5,' INEQUIV. ATOMS/UC USING ', &
       & A4,'. AUGMENTATION')") &
       LATTIC,NAT,RelMode
  WRITE(unit_out,"(/' LATTIC CONSTANTS:',3F12.7,' (in bohr)' &
       & /' UNIT CELL ANGLES:',3F12.7,' (in degree)')") &
       A,ALPHA,BETA,GAMMA
  CALL LATGEN(LATTIC,A,ALPHA,BETA,GAMMA,ORTHO,PRIM)

  NMT=0
  INDEX=0
  atoms: DO JATOM=1,NAT
     INDEX1 = INDEX+1
     INDEX=INDEX+1

1030 FORMAT(4X,I4,1X,3(3X,F10.7):/15X,I2)
     READ(unit_struct,1030) IATNR(INDEX),(POS(I,INDEX),I=1,3),MULT(JATOM)

     IATNR(INDEX) = ABS(IATNR(INDEX))

     DO MU=2,MULT(JATOM)
        INDEX=INDEX+1
        READ(unit_struct,1030) IATNR(INDEX),(POS(I,INDEX),I=1,3)
        IATNR(INDEX) = ABS(IATNR(INDEX))
     end DO

     READ(unit_struct, '(A10,5X,I5,2(5X,F10.5),5X,F5.2)') ANAME,JRI(JATOM),RNOT(JATOM),RMT(JATOM), &
          ZNUC(JATOM)
     DX(JATOM) = LOG( RMT(JATOM)/RNOT(JATOM) ) / DBLE( JRI(JATOM)-1 )
     IF(JRI(JATOM).GT.NRAD) call croak('NRAD too small')

     DO I=1,JRI(JATOM)
        RM(I,JATOM) = RNOT(JATOM) * EXP( DBLE(I-1)*DX(JATOM) )
     end DO

     RMT(JATOM) = RM(JRI(JATOM),JATOM)
     READ(unit_struct, '(20X,3F10.8)') ((ROTLOC(I,J,INDEX1),I=1,3),J=1,3)

     DO J=1,NMT
        IF(RMT(JATOM).EQ.RAD(J))THEN
           IRAD(JATOM) = J
           cycle atoms
        ENDIF
     end DO

     NMT = NMT + 1
     RAD(NMT) = RMT(JATOM)
     IRAD(JATOM) = NMT
  end DO atoms

  NDAT = INDEX
  allocate ( AUG(NRAD,(LMAX7+1)*(LMAX7+1),NDat))

  write(unit_out, &
       "(/' SYMMETRY OPERATIONS'                                     &
       & /' -------------------'                                     &
       & /' y = {Q|t}(x) : y_i = Sum(j) Q_ij x_j + t_i',             &
       & /' with y_i and x_j in conventional fractional coordinates' &
       &//' symm    Q(:,1)  Q(:,2)  Q(:,3)    t(:)')")
  READ(unit_struct, '(I4)') IORD
  IF(IORD.GT.NSYM) call croak('NSYM too small')

  DO K=1,IORD
     READ(unit_struct, '(3(3I2,F10.5/))') ((IMAT(J,I,K),J=1,3),TAU(I,K),I=1,3)
     write(unit_out,"(I5,3I8,F13.5/2(5X,3I8,F13.5/))") &
          K,((IMAT(J,I,K),J=1,3),TAU(I,K),I=1,3)
  end DO

  write(unit_out, &
       "( ' POSITIONS OF THE BASIS ATOMS'                &
       & /' ----------------------------'                &
       & /' x = Sum(j=1,3) f_i a_i  with  f_i in [0,1[', &
       &  '  (in the conventional unit cell)'            &
       &//' atom    f_1      f_2      f_3   ',           &
       &  '        x [bohr]     y [bohr]     z [bohr]')")

  DO INDEX=1,NDAT
     DO J=1,3
        R(J) = POS(1,INDEX)*BR1(1,J) + POS(2,INDEX)*BR1(2,J) &
             + POS(3,INDEX)*BR1(3,J)
     end DO
     write(unit_out,"(I5,1X,3F9.5,3X,3F13.7)") &
          IATNR(INDEX),(POS(I,INDEX),I=1,3),(R(J),J=1,3)
  end DO

  IF(.NOT.PRIM) CALL trans(NDAT,POS,IORD,IMAT,TAU)
  CALL ROTDEF(NAT,MULT,POS,IOP)
!
!     << set up local rotations for all basis atoms (in Cartesian coord.) >>
  INDEX=0
  DO JATOM=1,NAT
     INDEX1 = INDEX+1
!       << store ROTLOC(1.atom) before updating it >>
     ROT0(:,:) = ROTLOC(:,:, index1)

     DO MU=1,MULT(JATOM)
        INDEX = INDEX+1
        CALL LOCDEF(ROT0,IMAT(1,1,IOP(INDEX)),ROTLOC(1,1,INDEX))
     end DO
  end DO
  
  write(unit_out, &
       "(/' SYMMETRY ADAPTED LOCAL ROTATION MATRICES' &
       & /' ----------------------------------------' &
       & /' x'' = T^-1(x) : x''_i = Sum(j) x_j T_ji'  &
       & ,'  with x_j and x''_i in Cartesian coord.'  &
       &//' atom  symm.      T(:,1)      T(:,2)      T(:,3)')")

  DO INDEX=1,NDAT
     write(unit_out,"(2(I5,1X),3F12.6,2(/12X,3F12.6)/)") &
          IATNR(INDEX),IOP(INDEX), &
          ((ROTLOC(I,J,INDEX),I=1,3),J=1,3)

     IF(.NOT.ORTH(ROTLOC(1,1,INDEX)))THEN
        WRITE(unit_out,"(/'CURRENT ROTATION MATRIX IS NOT ORTHOGONAL')")
        call croak('MATRIX NOT ORTHOGONAL')
     ENDIF
  end DO

!     << read in the evaluation grid >>
! -------------------------------------------------------------------------
! NPG        -- the total number of grid points of the evaluation grid
! RGRID(:,i) -- the i-th grid point in primitive fractional coordinates
! -------------------------------------------------------------------------
  CALL GRDGEN(GridMode,NP)
  allocate (psi(npg))
  allocate (tot(npg))
  tot=0

!     << find the surrounding primitve cell SPC of each atomic sphere >>
! ---------------------------------------------------------------
! SPC := Sum(i=1,3) [-s_i,+s_i] * a_i  with  s_i > 0
! and {a_1,a_2,a_3} being the primitive lattice vectors 
! ---------------------------------------------------------------
  CALL SPCGEN(NAT,RMT,ATMS)
  write(unit_out, &
       "(/' PRIMITIVE CELLS SURROUNDING THE MUFFIN TIN SPHERES' &
       & /' --------------------------------------------------' &
       & /' cell = Sum(i=1,3) [-s_i,+s_i] a_i',                 &
       &  '  with  a_i = primitive lattice vector'              &
       &//' atom     RMT         s_1      s_2      s_3')")

  DO JATOM=1,NAT
     write(unit_out,"(I5,F10.5,3X,3F9.5)") &
          JATOM,RMT(JATOM),(ATMS(I,JATOM),I=1,3)
  end DO
!     << determine into which sphere each grid point r falls >>
! -------------------------------------------------------------------------
! if in interstitial:
! RGRID(1:3,i) -- the i-th grid point in (global) Cartesian coordinates
!
! if in muffin tin sphere around R(NX,NY,NZ) + R0(IAT)
! IREG (  i) -- the atom IAT for the i-th grid point
! ILAT (:,i) -- the lattice vector (NX,NY,NZ) for the i-th grid point
!               in primitive fractional coordinates 
! IRI  (  i) -- the radial interval the i-th grid point falls in
! RGRID(:,i) -- the value r - R - R0 for the i-th grid point 
!               in the local Cartesian coordinates of atom IAT
! -------------------------------------------------------------------------
  NREG(:) = 0

  DO IG=1,NPG
     CALL FINDMT(RGRID(1,IG),ATMS,nat,NDAT,INDEX,ILAT(1,IG),IRI(IG),R)
     IREG(IG) = INDEX
     NREG(INDEX)= NREG(INDEX) + 1
     IF(INDEX.GT.0)THEN
!         << in muffin tin sphere around R(NX,NY,NZ) + R0(IAT) >>
!         << transform r-R-R0 into local Cartesian coordinates >>
! -----------------------------------------------------------------------
! psi(r) = psi'(T^-1(r)) with psi' as the (LAPW1)-provided wave functions
! -----------------------------------------------------------------------
!         << r' := T^-1(r) >>
        DO I=1,3
           RGRID(I,IG) = ROTLOC(I,1,INDEX)*R(1) &
                + ROTLOC(I,2,INDEX)*R(2) &
                + ROTLOC(I,3,INDEX)*R(3)
        end DO
     ELSE
!         << in interstitial >>
        RGRID(:,IG) = R(:)
     ENDIF
  end DO

  IF(NPG.GT.0)THEN
     write(unit_out, &
          "(/' GRID POINT DISTRIBUTION'               &
          & /' -----------------------'               &
          & /' region  atom  grid points  percentage' &
          & /I7,3X,'int',I10,2P,F11.1,' %')")         &
          0,NREG(0),NREG(0)/dble(NPG)

     DO INDEX=1,NDAT
        write(unit_out,"(I7,I5,I11,2P,F11.1,' %')") &
             INDEX,IATNR(INDEX),NREG(INDEX), &
             NREG(INDEX)/dble(NPG)
     end DO
  ENDIF

!     << finally transform POS into global Cartesian coord. >>
  DO INDEX=1,NDAT
     DO J=1,3
        R(J) = POS(1,INDEX)*BR2(1,J) + POS(2,INDEX)*BR2(2,J) &
             + POS(3,INDEX)*BR2(3,J)
     end DO

     POS(:,INDEX) = R(:)
  end DO

! #######################################
! # PART III : wave function processing #
! #######################################

!     << writing to output is done later >>
!
!     << set up the augmentation functions  >>
!     << for each group of equivalent atoms >>
  call paropen(unit_vector, vecfn, iproc, 1, FORM='unformatted')
  CALL AUGGEN(REL,NAT,WHPSI)
  close(unit_vector)
!
!     << write wave function options to output >>
  write(unit_out, &
       "(/' WAVE FUNCTION OPTIONS'                            &
       & /' ---------------------'                            &
       & /' evaluation of the ',A5,' relativistic component', &
       &  ' of the wave functions')")                         &
       WHPSI

  IF(DEPHAS) WRITE(unit_out,"( ' after dephasing of the wave function')")
  write(unit_out,"( ' data provided    :  ',A,'  [in ',A3,' units]')") &
       SWITCH//'(PSI)',IUNIT
  write(BUF(1:67), &
       "('k =..... , n =...... : ',A8,' [',A3,' units] ', &
       & '-- ',A5,' rel. component')")                    &
       SWITCH//'(PSI)',IUNIT,WHPSI

!     << set up constants and prefactors >>
! -----------------------------------------
!     prefactor = 1/sqrt( vol(UC) )
! -----------------------------------------
  PREFAC = PREFAC / SQRT(VUC)
!
!     << scale reciprocal basis vectors with 2pi >>
  BR3 = 2*PI*BR3
  BR4 = 2*PI*BR4

!     ---------
!     MAIN LOOP
!     ---------
  if(GridMode(1:3) == 'ANY') then
     write(unit_out,2150)'according to the grid point input file'
  elseif(GridMode(1:2) == '2D') then
     write(unit_out,2150) '((psi(ix,iy),ix=1,nx),iy=1,ny)'
  elseif(GridMode(1:2) == '3D') then
     write(unit_out,2150) '(((psi(ix,iy,iz),ix=1,nx),iy=1,ny),iz=1,nz)'

2150 format(/' ==================' &
          & /' WAVE FUNCTION DATA' &
          & /' ==================' &
          & /' order: ',A/)
  endif

!     << read in k points >>
  KKK = 0
  nmat = 0

  vecfiles: do iloop = 1, max(iproc,1)
     call paropen(unit_vector, vecfn, iproc, iloop, FORM='unformatted')
     ! read off two records from the fresh ‘vector’ file
     do i=1,NAT
        read(unit_vector)
        read(unit_vector)
     end do
     kpoints: do
        READ(unit_vector, IOSTAT=IERR)&
             (FK(I),I=1,3),BNAME,NMATK,NE,WEIGHT,KKKABS
        IF(IERR.LT.0) exit kpoints

        KKK = KKK + 1

        do j=1, chk%num_bands
           local_contribution = dcmplx(0d0,0d0)
           if (chk%have_disentangled) then
              do j2=1, chk%num_wann
                 if (rotate) then
                    do j3=1, chk%num_wann
                       local_contribution = local_contribution +&
                            rotmatrix(iseig,j3) * &
                            chk%u_matrix(j2,j3,kkk)*chk%u_matrix_opt(j,j2,kkk)
                    enddo
                 else
                    local_contribution = local_contribution + &
                         chk%u_matrix(j2,ISEIG,kkk)*chk%u_matrix_opt(j,j2,kkk)
                 endif
              enddo
           else
              if (rotate) then
                 do j2=1, chk%num_wann
                    local_contribution = local_contribution+ &
                         rotmatrix(j2,iseig)*chk%u_matrix(j,j2,kkk)
                 enddo
              else
                 local_contribution = chk%u_matrix(j,iseig,kkk)
              endif
           endif
           uu(j)=local_contribution
        enddo
        
        IF(IERR.GT.0) KKKABS = KKK
        NPW = NMATK - NLO

        if (nmatk > nmat) then
           nmat = nmatk

!           if (allocated(kx)) deallocate(kx,ky,kz, fj,dfj, bk,coef,cf)
           if (allocated(kx)) deallocate(kx,ky,kz, bk,coef,cf)

           allocate (kx(nmat),ky(nmat),kz(nmat))
!           allocate ( FJ(0:LMAX7,NMAT,NMT),DFJ(0:LMAX7,NMAT,NMT))
           allocate ( BK(3,NMAT),COEF(NMAT),CF(NMAT))
        end if
        !! for some reason, different results if [D]FJ not
        !! reallocated to the new NMATK every time
        if(allocated(fj)) deallocate(fj,dfj)
        allocate ( FJ(0:LMAX7,NMATk,NMT),DFJ(0:LMAX7,NMATk,NMT))

        !       << read in PW and local orbital wave vectors K >>
        READ(unit_vector) (KX(IK),KY(IK),KZ(IK),IK=1,NMATK)
        !
        write(unit_out, "(' k-point',I4,' : ',A10,&
             &' = (',F8.4,',',F8.4,',',F8.4,')'/' ',55('-'))")&
             KKKABS,BNAME,(FK(I),I=1,3)


        !       << transform K+k into global Cartesian coordinates >>
        DO IK=1,NMATK
           BK(:,ik) = BR3(1,:)*(KX(ik)+FK(1)) + BR3(2,:)*(KY(ik)+FK(2)) + BR3(3,:)*(KZ(ik)+FK(3))
        end DO

        !       << transform k into primitive fractional coordinates >>
        IF(.NOT.PRIM)THEN
           DO J=1,3
              R(J) = BR3(1,J)*FK(1) + BR3(2,J)*FK(2) + BR3(3,J)*FK(3)
           end DO
           DO I=1,3
              FK(I) = BR2(I,1)*R(1) + BR2(I,2)*R(2) + BR2(I,3)*R(3)
           end DO
        ELSE
           FK(:) = 2*PI * FK(:)
        ENDIF

!       << load spherical Bessel functions and their derivatives >>
! -----------------------------------------------------------------
! FJ (l,:,Rmt) :      j_l(|K+k|*r) at r = Rmt_a for all PW's K+k
! DFJ(l,:,Rmt) : d/dr j_l(|K+k|*r) at r = Rmt_a for all PW's K+k
! -----------------------------------------------------------------
        call BESSEL(NPW,NMATK,BK,NMT,RAD,LMAX7,FJ,DFJ)
!
!       << read in individual eigen states >>
       
        COEF=0
        DO I=1,NE
340        READ(unit_vector,IOSTAT=IERR) IE,E,IEABS
           READ(unit_vector) (CF(IK),IK=1,NMATK)
           IF (ie>=nemin .and. ie<=nemax) THEN
              do ik=1,nmatk
                 COEF(ik)=COEF(ik)+uu(ie-nemin+1)*cf(ik)
              enddo
           ENDIF
        ENDDO
!
!         << normalize the eigen functions properly >>
        COEF = PREFAC*COEF
!
!         << set up the eigen state's augmentation coefficients >>
!         << ALM, BLM, and CLM for each atom in the unit cell   >>
! ----------------------------------------------------------------------------
! in the muffin thin sphere around R + R_a one has:
!
! psi(r) = e^ikR Sum(lm) w_lm,a(|r-R-R_a|) Y*(*)_lm(T_a^-1(r-R-R_a))   with
!
! w_lm,a(r) = 4pi*i^l [ A_l,m,a *      u_l,a(r,E _l,a) +
!                       B_l,m,a * d/dE u_l,a(r,E _l,a) +
!                       C_l,m,a *      u_l,a(r,E'_l,a) ] * Rmt_a^2
!
! Here (*) stands for an optional additional complex conjugation on Y*_lm(...)
! WIEN95 : (*) =     and hence   *(*) = *
! WIEN97 : (*) = *   and hence   *(*) = 
! ----------------------------------------------------------------------------
        IL1 = NPW
        LATOM  = 0
        DO JATOM=1,NAT
           IMAX = JRI(JATOM)
           DO JNEQ=1,MULT(JATOM)
              LATOM = LATOM + 1
              CALL AUGPW(LATOM,NPW,ALM,ROTLOC,Y,bk,coef,nmatk)

              IL=IL1
              CALL AUGLO(LATOM,IL,ALM,ROTLOC,Y,bk,coef,nmatk)
              IF (JNEQ.EQ.MULT(JATOM)) IL1=IL
                          
!             << add 4pi*i^l Rmt^2 factor to ALM, BLM and CLM and >>
!             << set up the eigen state's augmentation functions  >>
! ----------------------------------------------------------------------------
! w_lm,a(r) = 4pi*i^l Rmt_a^2 * [ A_l,m,a *      u_l,a(r,E _l,a) + 
!                                 B_l,m,a * d/dE u_l,a(r,E _l,a) +
!                                 C_l,m,a *      u_l,a(r,E'_l,a) ]
! ----------------------------------------------------------------------------
              PHS = 4*PI * RMT(JATOM)*RMT(JATOM)
              LM = 0
              DO L=0,LMAX7
                 DO M=-L,L
                    LM=LM+1
                    DO I=1,IMAX
                       ALMSUM=(0.d0,0.d0)
                       DO irf=1,nrf
                          ALMSUM=ALMSUM+ALM(LM,irf)*RRAD(I,L,irf,JATOM)
                       ENDDO
                       AUG(I,LM,LATOM)=ALMSUM*PHS
                    end DO
                 end DO
                 PHS = (0.0D0,1.0D0)*PHS
              end DO
           end DO
        end DO

!
!         << now evaluate the wave function on the grid >>
        DO IG=1,NPG
           IF(IREG(IG).EQ.0)THEN
!             << grid point in interstitial >>
              IF(LARGE)THEN
                 CALL WAVINT(RGRID(1,IG),NPW,PSI(IG),bk,coef,nmatk)
              ELSE
                 PSI(IG)=(0.0D0,0.0D0)
              ENDIF
           ELSE
!             << grid point in atomic sphere at R(ILAT) + R_a(IREG) >>
!             << load calc. Bloch factor e^ikR >>
              ARG = FK(1)*ILAT(1,IG) + FK(2)*ILAT(2,IG)  &
                   + FK(3)*ILAT(3,IG)
              PHS = DCMPLX(COS(ARG),SIN(ARG))
              CALL WAVSPH(RGRID(1,IG),PHS,IREG(IG),IRI(IG),PSI(IG),Y)
           ENDIF
!          write(93,*)IG,IREG(IG),PSI(IG)
        end DO

        DO IG=1,NPG
           tot(ig)=tot(ig)+psi(ig)
        ENDDO
        cycle kpoints

!         << correct for averaged phase factor >>
        IF(DEPHAS)THEN
           PHSAV = 0.0D0
           DO IG=1,NPG
              IF(ABS(PSI(IG)).GT.1D-18) &
                   PHSAV =  &
                   PHSAV + DMOD(DATAN2(DIMAG(PSI(IG)),DBLE(PSI(IG)))+PI,PI)
           end DO
           PHSAV = PHSAV/DBLE(NPG)
           PHS   = DCMPLX(COS(PHSAV),-SIN(PHSAV))
           DO IG=1,NPG
              PSI(IG) = PSI(IG)*PHS
           end DO
        ENDIF

        WRITE(BUF(1:20),"('k =',I5,' , n =',I6)") KKKABS,IEABS
        WRITE(unit_psink,'(A)') BUF(1:67)

!         << echo wave functions to the output >>
        write(unit_out, &
             "( ' band',I4,' : E = ',F12.6,'  (k-point:',I4,')')") &
             IEABS,E,KKKABS

        if(IE.lt.NE) goto 340
     end do kpoints
  end do vecfiles

  write(unit_psink, 3020)(abs(tot(IG))**2/kkk,IG=1,NPG)
  write(unit_psiarg,3020)(atan2(dimag(tot(IG)),real(tot(IG))),IG=1,NPG)
3020 format(1P,10E16.8)

  RETURN


contains
  subroutine wfopen(unit, file, iwf, iostat, status, form)
    integer,      intent(in)            :: unit
    character(*), intent(in)            :: file
    integer,      intent(in)            :: iwf
    integer,      intent(out), optional :: iostat
    character(*), intent(in),  optional :: status
    character(*), intent(in),  optional :: form

    character(len=BUFSZ) :: tmpnam
    character(len=   11) :: f = 'FORMATTED'
    character(len=    7) :: s = 'UNKNOWN'
    integer              :: idot

    if (present(status)) s = status
    if (present(form))   f = form

    if (iwf==0) then
       tmpnam = file
    else
       idot = scan(file, '.', BACK=.true.)

       if (idot/=0) then
          write(tmpnam, '(A, "_", I0, A)') &
               file(1:idot-1), iwf, trim(file(idot:))
       else
          write(tmpnam, '(A, "_", I0)') trim(file), iwf
       end if
    end if

    if (present(iostat)) then
       open(unit, FILE=tmpnam, STATUS=s, FORM=f, IOSTAT=iostat)
    else
       open(unit, FILE=tmpnam, STATUS=s, FORM=f)
    end if
  end subroutine wfopen
end subroutine main


!!/---
!! Local Variables:
!! mode: f90
!! End:
!!\---
!!
!! Time-stamp: <2015-12-22 20:16:28 assman@faepop36.tu-graz.ac.at>
