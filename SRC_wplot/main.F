!!! wien2wannier/SRC_wplot/main.F

subroutine main()
  use struct,    only: pos, RMT
  use radgrd,    only: rm, dx
  use latt,      only: br1, br2, br3, br4
  use sym2,      only: tau, imat, iord
  use lolog,     only: lapw, ilo, Nlo
  use loabc,     only: Alo
  use atspdt,    only: p, dp
  use radfu,     only: rrad
  use bessfu,    only: fj, dfj, rad, irad
  use work,      only: aug
  use grid,      only: rgrid, ireg, ilat, iri, npg
  use const,     only: DPk, BUFSZ, PI, error_unit
  use util,      only: paropen, uppercase
  use clio,      only: croak, carp
  use Wannier90, only: chk_t, chk_read
  use structmod, only: struct_t, struct_read
  use param,     only: Nrad, Nsym, NRF, NLOAt, LOmax, Lmax7, idx_wann,      &
       unit_in, unit_out, unit_struct, unit_vector, unit_chk, unit_rot,     &
       unit_inwf, unit_psink, unit_psiarg, outfn, vecfn, psinkfn, psiargfn, &
       iproc

  implicit none

  interface
     logical function orth(A)
       use const, only: DPk

       real(DPk), intent(in) :: A(3,3)
     end function orth
  end interface

  character(BUFSZ) :: buf
  character( 3)    :: switch, iunit
  character( 4)    :: GridMode
  character( 5)    :: whpsi
  character(10)    :: bname

  logical      :: LARGE
  complex(DPk) :: ALMSUM

  type(struct_t) :: stru
  type(chk_t)    :: chk

#ifndef _COMPLEX_
  REAL(DPk),    allocatable :: CF(:)
#else
  COMPLEX(DPk), allocatable :: CF(:)
#endif

  COMPLEX(DPk), allocatable ::  COEF(:) !changed by pwissgott

  real(DPk),allocatable :: BK(:,:)

  COMPLEX(DPk),allocatable :: PSI(:),tot(:)
  COMPLEX(DPk) ALM((LMAX7+1)*(LMAX7+1),NRF), Y((LMAX7+1)*(LMAX7+1))
  COMPLEX(DPk) PHS
  real(DPk),allocatable ::  ROTLOC(:,:,:),atms(:,:)
  integer,allocatable :: IOP(:),nreg(:)
  integer,allocatable ::  KX(:),KY(:),KZ(:)
  complex(DPk), allocatable :: uu(:), rotmat(:,:)
  real(DPk),allocatable :: tmprot(:)
  LOGICAL :: DEPHAS, rotate=.false.

  integer :: i, ie, ieabs, ierr, iG, ik, il, il1, imax, index, index1, irf
  integer :: iseig, iloop
  integer :: j, j2, j3, jneq, jmult, k, kkk, kkkabs, l, latom, lm, m, mu
  integer :: Nb, NE, NEmax, NEmin, Nmatk, Nmat, Nmt, Npw
  integer :: NP(3), Nskip

  real(DPK) :: ROT0(3,3), R(3), FK(3), phsav, prefac
  real(DPK) :: arg, E, weight

  complex(DPk) :: local_contribution

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
  read (unit_in,*) iseig, rotate
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
     write(unit_out,*) "rotate with matrix from unit ", unit_rot
     allocate(rotmat(chk%num_wann, chk%num_wann), tmprot(chk%num_wann*2))

     do i=1,chk%num_wann
        read(unit_rot,*) tmprot

        do j=1,chk%num_wann
           rotmat(i,j) = cmplx(tmprot(2*j-1), tmprot(2*j), DPk)
        enddo
     enddo
     close(unit_rot)

     if (any(abs( matmul(rotmat, conjg(transpose(rotmat))) &
          &       - Id(chk%num_wann))                      &
          &  > 1e-10)) then
        call carp("rotation matrix does not look unitary")
        write(buf, '( "(", I0, "(F10.5, '' + '', F7.5, ''i''))" )') &
             chk%num_wann
        write(error_unit, buf) transpose(rotmat)
     end if
  endif

  rewind(unit_in)               ! for grdgen()
  close(unit_inwf)

! ###################################
! # PART II : structural processing #
! ###################################
!
!     << set up the lattice >>
  call struct_read(unit_struct, stru)
  close(unit_struct)

  allocate ( POS(3,stru%Nneq*48),             RMT(stru%Nneq), DX(stru%Nneq))
  allocate ( ILO(0:LOMAX, stru%Nneq),         LAPW(0:LMAX7,stru%Nneq)      )
  allocate ( ALO(0:LOMAX,NLOAT,NRF,stru%Nneq),RM(NRAD,stru%Nneq)           )
  allocate ( P(0:LMAX7,NRF,stru%Nneq),        DP(0:LMAX7,NRF,stru%Nneq)    )
  allocate ( RRAD(NRAD,0:LMAX7,NRF,stru%Nneq)                              )
  allocate ( RAD(stru%Nneq),                  IRAD(stru%Nneq)              )
  allocate ( ROTLOC(3,3,stru%Nneq*48),        atms(3,stru%Nneq)            )
  allocate ( IOP(stru%Nneq*48),               nreg(0:stru%Nneq*48)         )

  write(unit_out, '(/A/)') trim(stru%title)
  write(unit_out,"(1X,A4,' LATTIC WITH ',I5,' INEQUIV. ATOMS/UC USING ', &
       & A4,'. AUGMENTATION')") &
       stru%lattic, stru%Nneq, stru%mode
  write(unit_out,"(/' LATTIC CONSTANTS:',3F12.7,' (in bohr)' &
       & /' UNIT CELL ANGLES:',3F12.7,' (in degree)')") &
       stru%a, stru%alpha(:)

  call latgen(stru)

  rmt  = stru%rmt
  dx   = log(stru%rmt / stru%r0) / real(stru%npt-1, DPk)

  ! FIXME: Nrad could be replaced by max(stru%npt)
  if (any(stru%Npt > Nrad)) call croak('NRAD too small')

  NMT = 0
  atoms: do jneq=1,stru%Nneq
     do i=1, stru%Npt(jneq)
        RM(i, jneq) = stru%r0(jneq) * exp( real(i-1, DPk) * dx(jneq) )
     end do

     RMT(jneq) = RM(stru%Npt(jneq), jneq)


     rotloc(:,:, sum(stru%mult(1:jneq-1))+1) = &
          transpose(stru%locrot(:,:, jneq))

     do j=1,NMT
        if(RMT(jneq) == rad(j))then
           irad(jneq) = j
           cycle atoms
        ENDIF
     end DO

     NMT        = NMT + 1
     rad(NMT)   = RMT(jneq)
     irad(jneq) = NMT
  end do atoms

  allocate ( AUG(NRAD,(LMAX7+1)*(LMAX7+1), stru%nat))

  write(unit_out, &
       "(/' SYMMETRY OPERATIONS'                                     &
       & /' -------------------'                                     &
       & /' y = {Q|t}(x) : y_i = Sum(j) Q_ij x_j + t_i',             &
       & /' with y_i and x_j in conventional fractional coordinates' &
       &//' symm    Q(:,1)  Q(:,2)  Q(:,3)    t(:)')")

  iord = size(stru%rsym, 3)
  IF(IORD.GT.NSYM) call croak('NSYM too small')

  imat = stru%rsym
  tau  = stru%rtrans
  DO K=1,IORD
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

  do INDEX=1,stru%nat
     do J=1,3
        R(J) = stru%pos(1,index)*BR1(1,J) + stru%pos(2,index)*BR1(2,J) &
             + stru%pos(3,index)*BR1(3,J)
     end do
     write(unit_out,"(I5,1X,3F9.5,3X,3F13.7)") &
          stru%Neq2at(INDEX),(stru%pos(I,INDEX),I=1,3),(R(J),J=1,3)
  end do

  ! trans() updates stru%pos, IMAT, TAU
  if (.not. stru%prim) call trans(stru)

  ! rotdef() updates POS, sets IOP
  call rotdef(stru, IOP)

!     << set up local rotations for all basis atoms (in Cartesian coord.) >>
  INDEX=0
  DO JNEQ=1,stru%Nneq
     INDEX1 = INDEX+1
!       << store ROTLOC(1.atom) before updating it >>
     ROT0(:,:) = ROTLOC(:,:, index1)

     DO MU=1,STRU%MULT(JNEQ)
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

  DO INDEX=1,stru%nat
     write(unit_out,"(2(I5,1X),3F12.6,2(/12X,3F12.6)/)") &
          STRU%NEQ2AT(INDEX),IOP(INDEX), &
          ((ROTLOC(I,J,INDEX),I=1,3),J=1,3)

     IF(.NOT.ORTH(ROTLOC(:,:,INDEX)))THEN
        WRITE(unit_out,"(/'CURRENT ROTATION MATRIX IS NOT ORTHOGONAL')")
        call croak('Local rotation matrix not orthogonal')
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
  CALL SPCGEN(stru%Nneq,RMT,ATMS)
  write(unit_out, &
       "(/' PRIMITIVE CELLS SURROUNDING THE MUFFIN TIN SPHERES' &
       & /' --------------------------------------------------' &
       & /' cell = Sum(i=1,3) [-s_i,+s_i] a_i',                 &
       &  '  with  a_i = primitive lattice vector'              &
       &//' atom     RMT         s_1      s_2      s_3')")

  DO JNEQ=1, stru%Nneq
     write(unit_out,"(I5,F10.5,3X,3F9.5)") &
          JNEQ,RMT(JNEQ),(ATMS(I,JNEQ),I=1,3)
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
     call findmt(RGRID(1,IG), atms, stru, index, iLat(1,ig), iri(ig), R)
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
          0,NREG(0),NREG(0)/real(NPG, DPk)

     DO INDEX=1,stru%nat
        write(unit_out,"(I7,I5,I11,2P,F11.1,' %')") &
             INDEX,STRU%NEQ2AT(INDEX),NREG(INDEX), &
             NREG(INDEX)/real(NPG, DPk)
     end DO
  ENDIF

!     << finally transform POS into global Cartesian coord. >>
  DO INDEX=1,stru%nat
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
  call auggen(stru%mode=='RELA', stru, whpsi)
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
  PREFAC = PREFAC / SQRT(STRU%VOL)
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
     do i=1,stru%Nneq
        read(unit_vector)
        read(unit_vector)
     end do
     kpoints: do
        ! FIXME: every one of these reads fails with “input statement
        ! requires too much data
        read(unit_vector, IOSTAT=ierr) FK, BNAME, NMATK, NE, WEIGHT, KKKABS
        if (ierr<0) exit kpoints

        KKK = KKK + 1

        do j=1, chk%num_bands
           local_contribution = 0
           if (chk%have_disentangled) then
              do j2=1, chk%num_wann
                 if (rotate) then
                    do j3=1, chk%num_wann
                       local_contribution = local_contribution +&
                            rotmat(iseig,j3) * &
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
                         rotmat(j2,iseig)*chk%u_matrix(j,j2,kkk)
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
        if(.not. stru%prim)then
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
        do I=1,NE
           read(unit_vector,IOSTAT=IERR) IE,E,IEABS
           read(unit_vector) (CF(IK),IK=1,NMATK)
           if (ie>=nemin .and. ie<=nemax) then
              do ik=1,nmatk
                 COEF(ik)=COEF(ik)+uu(ie-nemin+1)*cf(ik)
              end do
           end if
        end do
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
        DO JNEQ=1,stru%Nneq
           IMAX = STRU%NPT(JNEQ)
           DO JMULT=1,STRU%MULT(JNEQ)
              LATOM = LATOM + 1
              call augpw(latom, Npw, Alm, rotloc, Y, bk, coef, Nmatk, &
                   stru%neq2at)

              IL=IL1
              call auglo(latom, il, Alm, rotloc, Y, bk, coef, nmatk, &
                   stru)
              IF (JMULT.EQ.STRU%MULT(JNEQ)) IL1=IL

!             << add 4pi*i^l Rmt^2 factor to ALM, BLM and CLM and >>
!             << set up the eigen state's augmentation functions  >>
! ----------------------------------------------------------------------------
! w_lm,a(r) = 4pi*i^l Rmt_a^2 * [ A_l,m,a *      u_l,a(r,E _l,a) +
!                                 B_l,m,a * d/dE u_l,a(r,E _l,a) +
!                                 C_l,m,a *      u_l,a(r,E'_l,a) ]
! ----------------------------------------------------------------------------
              PHS = 4*PI * RMT(JNEQ)*RMT(JNEQ)
              LM = 0
              DO L=0,LMAX7
                 DO M=-L,L
                    LM=LM+1
                    DO I=1,IMAX
                       ALMSUM=0
                       DO irf=1,nrf
                          ALMSUM=ALMSUM+ALM(LM,irf)*RRAD(I,L,irf,JNEQ)
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
              if(LARGE)then
                 call WAVINT(RGRID(1,IG),NPW,PSI(IG),bk,coef,nmatk)
              else
                 PSI(IG)=(0.0D0,0.0D0)
              endif
           ELSE
!             << grid point in atomic sphere at R(ILAT) + R_a(IREG) >>
!             << load calc. Bloch factor e^ikR >>
              ARG = FK(1)*ILAT(1,IG) + FK(2)*ILAT(2,IG)  &
                   + FK(3)*ILAT(3,IG)
              PHS = cmplx(cos(arg), sin(arg), DPk)
              call WAVSPH(Rgrid(1,ig), PHS, ireg(ig), iri(ig), Psi(ig), Y, &
                   stru%neq2at)
           ENDIF
        end DO

        DO IG=1,NPG
           tot(ig)=tot(ig)+psi(ig)
        ENDDO
     end do kpoints
  end do vecfiles

  !         << correct for averaged phase factor >>
  if(DEPHAS) then
     PHSAV = 0
     do IG=1,NPG
        if(abs(tot(IG)) > 1e-18_DPk) &
             PHSAV =  &
             PHSAV + mod(atan2(aimag(tot(ig)), real(tot(ig)))+PI, PI)
     end do
     PHSAV = PHSAV/real(NPG, DPk)
     PHS   = cmplx(cos(PHSAV), -sin(PHSAV), DPk)

     write(0,*) PHS

     tot = tot * PHS
  end if

  write(unit_psink, 3020) abs(tot/kkk)**2
  write(unit_psiarg,3020) atan2(aimag(tot), real(tot))
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

  pure function Id(N)
    integer, intent(in) :: N
    real(DPk)           :: Id(N,N)

    integer :: i

    Id = 0
    forall(i=1:N) Id(i,i) = 1
  end function Id
end subroutine main


!!/---
!! Local Variables:
!! mode: f90
!! End:
!!\---
!!
!! Time-stamp: <2016-02-10 17:25:08 assman@faepop36.tu-graz.ac.at>
