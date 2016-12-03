!!! wien2wannier/SRC_wplot/modules_rc.F
!!!
!!!    Modules for wien2wannier.  This file contains modules that
!!!    depend on real/complex compilation.

!!/=== Defined here: =============================
!!
!! Procedure module main_m
!!
!!\===============================================

module     main_m
  use const, only: DPk

  implicit none
  private
  public :: main

  real(DPk), parameter :: MAT_TOL = 1e-10_DPk, PHASE_TOL = 1e-10_DPk

  interface orth
     module procedure orthogonal, unitary
  end interface orth
contains
subroutine main()
  use sym2,       only: init_sym2, rtrans=>trans, imat
  use lolog,      only: lapw, ilo, Nlo
  use loabc,      only: Alo
  use atspdt,     only: p, dp
  use radfu,      only: rrad
  use bessfu,     only: fj, dfj, rad, irad
  use work,       only: aug
  use grid,       only: ireg, ilat, iri
  use const,      only: DPk, BUFSZ, TAU, error_unit
  use util,       only: paropen, uppercase
  use clio,       only: croak, carp
  use Wannier90,  only: chk_t, chk_read
  use structmod,  only: struct_t, struct_read
  use inwfmod,    only: inwf_t, inwf_read
  use inwplotmod, only: inwplot_t, inwplot_read
  use param,      only: NRF, NLOAt, LOmax, Nrad, wien2wannier_version, &
       &                unit_in, unit_out, unit_struct, unit_vector
  use wplot,      only: unit_inwf, unit_psink, unit_psiarg, unit_chk, &
       &                unit_rot, outfn, vecfn, psinkfn, psiargfn, gridfn, &
       &                iproc, Lmax7, idx_wann

  !! procedure includes
  use bessel_m
  use wavsph_m
  use wavint_m
  use rotdef_m
  use locdef_m
  use spcgen_m
  use findmt_m
  use auggen_m
  use auglo_m
  use augpw_m
  use trans_m

  implicit none

  character(BUFSZ) :: buf
  character(10)    :: bname

  type(struct_t)   :: stru
  type(chk_t)      :: chk
  type(inwf_t)     :: inwf
  type(inwplot_t)  :: inwp

#ifndef _COMPLEX_
  real(DPk),    allocatable :: CF(:)
#else
  complex(DPk), allocatable :: CF(:)
#endif

  complex(DPk), allocatable :: coef(:), psi(:), tot(:), uu(:), rotmat(:,:)
  real(DPk),    allocatable :: BK(:,:), atms(:,:), tmprot(:)
  real(DPk),    allocatable :: rotloc_sym(:,:,:), pos_sym(:,:)
  integer,      allocatable :: iop(:), nreg(:), KX(:),KY(:),KZ(:)

  complex(DPk) :: Alm((LMAX7+1)**2, NRF), Y((LMAX7+1)**2)
  complex(DPk) :: phs, AlmSum, local_contribution

  real(DPK)    :: R(3), FK(3), phsav, prefac, arg, E, weight

  integer      :: i, ie, ieabs, ierr, iG, ik, il, il1, imax, index
  integer      :: j, j2, j3, jneq, jmult, k, kkk, kkkabs, l, latom, lm, m, mu
  integer      :: Nb, NE, Nmatk, Nmat, Nmt, Npw, Npt, iloop, irf

  character, parameter :: xyz(3) = (/ 'x', 'y', 'z' /)

  character(*), parameter :: &
       fmt_grid  = "('# ', A1, i7, 2x, f13.7, 2(2x, f10.5))"

  call struct_read(unit_struct, stru)
  close(unit_struct)

  call inwplot_read(unit_in, inwp, stru)
  close(unit_in)

!!! Now we know the WF index INWP%WFIDX and can open files accordingly …
  call wfopen(unit_psink,  psinkfn,  inwp%WFidx, STATUS='unknown')
  call wfopen(unit_psiarg, psiargfn, inwp%WFidx, STATUS='unknown')
  call wfopen(unit_out,    outfn,    inwp%WFidx, STATUS='unknown')

!!! … and begin writing to unit_out.
  write(unit_out, '("wplot ", A /)') wien2wannier_version

  write(unit_out,"( ' selected Wannier function:',I5)") inwp%WFidx

  write(unit_out, '(/A/)') trim(stru%title)
  write(unit_out,"(1X,A4,' LATTIC WITH ',I5,' INEQUIV. ATOMS/UC USING ', &
       & A4,'. AUGMENTATION')") &
       stru%lattic, stru%Nneq, stru%mode
  write(unit_out,"(/' LATTIC CONSTANTS:',3F12.7,' [Bohr]' &
       & /' UNIT CELL ANGLES:',3F12.7,' [°]')") &
       stru%a, stru%alpha(:)

  write(unit_out,"(/' REAL SPACE LATTICE VECTORS a1, a2, a3 (in bohr)' &
       &           /' -----------------------------------------------' &
       &    /' CONVENTIONAL UNIT CELL :'/(' a_',I1,'    = ',3F12.7))") &
       ( i, stru%conv_dir(:,i), i=1,3 )
  write(unit_out,"(/' PRIMITIVE UNIT CELL :'/(' a_',I1,'    = ',3F12.7))") &
       ( i, stru%prim_dir(:,i), i=1,3 )
  write(unit_out,"(/' RECIPROCAL LATTIC VECTORS b1, b2, b3 (in 1/bohr)' &
       &           /' ------------------------------------------------' &
       &     /' CONVENTIONAL UNIT CELL :'/(' b_',I1,'/τ  = ',3F12.7))") &
       ( i, stru%conv_rec(:,i)/TAU, i=1,3 )
  write(unit_out,"(/' PRIMITIVE UNIT CELL :'/(' b_',I1,'/τ  = ',3F12.7))") &
       ( i, stru%prim_rec(:,i)/TAU, i=1,3 )

  write(unit_out, &
       "(/' SYMMETRY OPERATIONS'                                     &
       & /' -------------------'                                     &
       & /' y = {Q|t}(x) : y_i = Sum(j) Q_ij x_j + t_i',             &
       & /' with y_i and x_j in conventional fractional coordinates' &
       &//' symm    Q(:,1)  Q(:,2)  Q(:,3)    t(:)')")
  do k = 1, size(stru%rsym, 3)
     write(unit_out,"(I5,3I8,F13.5/2(5X,3I8,F13.5/))") &
          k,( stru%rsym(:,i,k), stru%rtrans(i,k),i=1,3)
  end do

  write(unit_out, &
       "( ' POSITIONS OF THE BASIS ATOMS'                &
       & /' ----------------------------'                &
       & /' x = Sum(j=1,3) f_i a_i  with  f_i in [0,1[', &
       &  '  (in the conventional unit cell)'            &
       &//' atom    f_1      f_2      f_3   ',           &
       &  '        x [Bohr]     y [Bohr]     z [Bohr]')")
  do INDEX=1,stru%nat
     write(unit_out,"(I5,1X,3F9.5,3X,3F13.7)") stru%Neq2at(index), &
          stru%pos(:,index), matmul(stru%conv_dir, stru%pos(:,index))
  end do

  Npt = size(inwp%rpts, 2)

  allocate( P   (      0:Lmax7, Nrf, stru%Nneq), rotloc_sym(3, 3, stru%Nat) )
  allocate( DP  (      0:Lmax7, Nrf, stru%Nneq),    pos_sym(3,    stru%Nat) )
  allocate( LAPW(      0:Lmax7,      stru%Nneq),  RAD (   stru%Nneq)        )
  allocate( Rrad(Nrad, 0:Lmax7, Nrf, stru%Nneq),  IRAD(   stru%Nneq)        )
  allocate( aug (Nrad,  (Lmax7+1)**2, stru%Nat),  IOP (   stru%Nneq*48)     )
  allocate( ALO (0:LOmax, NLOat, Nrf, stru%Nneq), nreg( 0:stru%Nneq*48)     )
  allocate( ILO (0:LOmax,             stru%Nneq), atms(3, stru%Nneq)        )
  allocate( psi(Npt), tot(Npt), iLat(3, Npt), iRi(Npt), iReg(Npt)           )
  tot=0

  pos_sym = stru%pos

  call init_sym2(stru)

  ! FIXME: Nrad could be replaced by max(stru%npt)
  if (any(stru%Npt > Nrad)) call croak('NRAD too small')

  ! trans() updates pos_sym, IMAT, RTRANS
  if (.not. stru%prim) call trans(pos_sym, stru)

  ! rotdef() updates pos_sym, sets IOP
  call rotdef(stru, iop, pos_sym)

  ! set up local rotations for all basis atoms (in Cartesian coord.)
  index=0
  do jneq=1,stru%Nneq
     do mu=1,stru%mult(jneq)
        index = index+1
        call locdef(stru%locrot(:,:,jneq), imat(:,:,iop(index)), &
             &      rotloc_sym (:,:,index), stru)
     end do
  end do

  write(unit_out, &
       "(/' SYMMETRY ADAPTED LOCAL ROTATION MATRICES' &
       & /' ----------------------------------------' &
       & /' x'' = T^-1(x) : x''_i = Sum(j) x_j T_ji'  &
       & ,'  with x_j and x''_i in Cartesian coord.'  &
       &//' atom  symm.      T(:,1)      T(:,2)      T(:,3)')")
  do INDEX=1,stru%Nat
     write(unit_out,"(2(I5,1X),3F12.6,2(/12X,3F12.6)/)") &
          stru%neq2at(index), iop(index), rotloc_sym(:,:,index)

     if (.not. orth(rotloc_sym(:,:,index))) then
        write(unit_out,"(/'CURRENT ROTATION MATRIX IS NOT ORTHOGONAL')")
        call croak('Local rotation matrix not orthogonal')
     endif
  end do

  extgrid: if (inwp%extgrid) then
     write(unit_psink,                                  &
          "('ANY   NP =',I0,/,                          &
          & 'order according to the grid points ',/,    &
          & 'provided in the input file ',/,            &
          & '`', A, '''')") &
          Npt, trim(gridfn)

     write(unit_out, &
          & "(' ARBITRARY LIST OF GRID POINTS' / &
          & ' -----------------------------' / &
          & ' number of grid points: ',I7)") Npt
  else
     write(unit_out, &
          &"(' 3D-NET OF GRID POINTS'                                    /&
          &' ---------------------'                                      /&
          &' number of grid points for x, y, z: ',3(I0,' '),              &
          &' (total: ',I0,')')")                                          &
          inwp%grid%N, Npt

     write(unit_out, "(/' PLOTTING AREA'                         &
          &            /' -------------'                         &
          &            /' x = Sum(j=1,3) f_i a_i  with  f_i in', &
          &             ' conventional fractional coordinates'   &
          &           //'           f_1      f_2      f_3   ',   &
          &             '        x [Bohr]     y [Bohr]     z [Bohr]')")

     write(unit_out, "(' origin', 1X, 3F9.5, 3X, 3F13.7)") &
          inwp%grid%frac(:,0), inwp%grid%cart(:,0)

     do i=1,3
        write(unit_out, "(' ',A1,'-axis',1X,3F9.5,3X,3F13.7)") &
             xyz(i), inwp%grid%frac(:,i), inwp%grid%cart(:,i)
     end do

     write(unit_out, "(/'         length [Bohr]')")
     do i=1,3
        write(unit_out, "(' ‖',A1,'‖ ',F13.7)") xyz(i), inwp%grid%len(i)
     end do

     write(unit_out, "(/'         cos(φ)        φ [°]')")
     j=3
     do i=1,2
        do j=i+1,3
           write(unit_out, "(' ∠(',A1,',',A1,') ',F10.7,3X,F8.3)") &
                xyz(i), xyz(j), cos(inwp%grid%ang(3)), inwp%grid%ang(3)
        end do
     end do

     write(unit_psink, "('# α     NP(α)   length(α)    ∠(α, x)    ∠(α, y)')")

     write(unit_psink, fmt_grid) &
          xyz(1), inwp%grid%N(1), inwp%grid%len(1)
     write(unit_psink, fmt_grid) &
          xyz(2), inwp%grid%N(2), inwp%grid%len(2), inwp%grid%ang(3)
     write(unit_psink, fmt_grid) &
          xyz(3), inwp%grid%N(2), inwp%grid%len(3), inwp%grid%ang(2:1:-1)

     write(unit_psink, &
          "('#order: (((w(x,y,z), z=1,#z), y=1,#y), x=1,#x)')")
  endif extgrid

  call inwf_read(unit_inwf, inwf)
  close(unit_inwf)

  Nb = inwf%bmax - inwf%bmin + 1

  allocate(uu(Nb))
  if (idx_wann.gt.0) then
     inwp%WFidx = idx_wann
  endif

!!! Read chk file for num_wann, u_matrix, and friends
  call chk_read(unit_chk, chk)
  close(unit_chk)

  write(unit_out,'("Read Wannier90 ‘chk’ file ", A)') trim(chk%header)

  if (chk%num_bands /= nb) &
       call croak("Error in read-in of chk file: &
       &inconsistency of number of Bloch bands")

  if (inwp%WFidx > chk%num_wann) &
       call croak("Error in read-in of chk file: &
          &number of Wannier functions is smaller than given index")

!!! Apply WF rotation
  if (inwp%WFrot) then
     write(unit_out,*) "rotate with matrix from unit ", unit_rot
     allocate(rotmat(chk%num_wann, chk%num_wann), tmprot(chk%num_wann*2))

     do i=1,chk%num_wann
        read(unit_rot,*) tmprot

        do j=1,chk%num_wann
           rotmat(i,j) = cmplx(tmprot(2*j-1), tmprot(2*j), DPk)
        enddo
     enddo
     close(unit_rot)

     if (.not. orth(rotmat)) then
        call carp("rotation matrix is not unitary")
        write(buf, '( "(", I0, "(F10.5, '' + '', F7.5, ''i''))" )') &
             chk%num_wann
        write(error_unit, buf) transpose(rotmat)
     end if
  endif

  NMT = 0
  atoms: do jneq=1,stru%Nneq
     do j=1,NMT
        if(stru%RMT(jneq) == rad(j))then
           irad(jneq) = j
           cycle atoms
        endif
     end do

     NMT        = NMT + 1
     rad(NMT)   = stru%RMT(jneq)
     irad(jneq) = NMT
  end do atoms

! find the surrounding primitve cell SPC of each atomic sphere
! ---------------------------------------------------------------
! SPC := Sum(i=1,3) [-s_i,+s_i] * a_i  with  s_i > 0
! and {a_1,a_2,a_3} being the primitive lattice vectors
! ---------------------------------------------------------------
  call SPCGEN(stru, atms)
  write(unit_out, &
       "(/' PRIMITIVE CELLS SURROUNDING THE MUFFIN TIN SPHERES' &
       & /' --------------------------------------------------' &
       & /' cell = Sum(i=1,3) [-s_i,+s_i] a_i',                 &
       &  '  with  a_i = primitive lattice vector'              &
       &//' atom     RMT         s_1      s_2      s_3')")

  do jneq=1, stru%Nneq
     write(unit_out,"(I5,F10.5,3X,3F9.5)") jneq, stru%RMT(jneq), atms(:,jneq)
  end do
! determine into which sphere each grid point r falls
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

  do ig=1,Npt
     call findmt(inwp%rpts(:,ig), atms, stru, pos_sym, index, &
          &      iLat(:,ig), iRi(ig), R)
     IREG(IG) = INDEX
     NREG(INDEX)= NREG(INDEX) + 1
     if(INDEX.gt.0)then
        ! in muffin tin sphere around R(NX,NY,NZ) + R0(IAT) transform
        ! r-R-R0 into local Cartesian coordinates psi(r) =
        ! psi'(T^-1(r)) with psi' as the (LAPW1)-provided wave
        ! functions

        ! r' := T^-1(r)
        inwp%rpts(:,ig) = matmul(rotloc_sym(:,:,index), r)
     else
        ! in interstitial
        inwp%rpts(:,ig) = r
     endif
  end do

  if(Npt.gt.0)then
     write(unit_out, &
          "(/' GRID POINT DISTRIBUTION'               &
          & /' -----------------------'               &
          & /' region  atom  grid points  percentage' &
          & /I7,3X,'int',I10,2P,F11.1,' %')")         &
          0,NREG(0),NREG(0)/real(Npt, DPk)

     do INDEX=1,stru%nat
        write(unit_out,"(I7,I5,I11,2P,F11.1,' %')") &
             INDEX,STRU%NEQ2AT(INDEX),NREG(INDEX), &
             NREG(INDEX)/real(Npt, DPk)
     end do
  endif

! finally transform POS into global Cartesian coord.
  do index=1,stru%Nat
     pos_sym(:,index) = matmul(stru%prim_dir, pos_sym(:, index))
  end do

! #######################################
! # PART III : wave function processing #
! #######################################

! set up the augmentation functions
! for each group of equivalent atoms
  call paropen(unit_vector, vecfn, iproc, 1, FORM='unformatted')
  call auggen(stru, inwp%large)
  close(unit_vector)

! write wave function options to output
  write(unit_out, &
       "(/' WANNIER FUNCTION OPTIONS'                         &
       & /' ------------------------'                         &
       & /' evaluation of the ',A,' relativistic component',  &
       &  ' of the wave functions in ',A,' units')")          &
       merge('large', 'small', inwp%large), trim(inwp%unit_name)

  if(inwp%dephas) &
       write(unit_out,"( ' after dephasing of the wave function')")

! set up constants and prefactors
  prefac = inwp%unit / sqrt(stru%vol)

!!! Main loop
  if(inwp%extgrid) then
     write(unit_out,2150)'according to the grid point input file'
  else
     write(unit_out,2150) '(((w(x,y,z), z=1,#z), y=1,#y), x=1,#x)'

2150 format(/' =====================' &
          & /' WANNIER FUNCTION DATA' &
          & /' =====================' &
          & /' order: ',A/)
  endif

! read in k points
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
                 if (inwp%WFrot) then
                    do j3=1, chk%num_wann
                       local_contribution = local_contribution &
                            + rotmat(inwp%WFidx,j3)         &
                            & * chk%u_matrix(j2,j3,kkk)        &
                            & * chk%u_matrix_opt(j,j2,kkk)
                    enddo
                 else
                    local_contribution = local_contribution  &
                         + chk%u_matrix(j2,inwp%WFidx,kkk) &
                         & * chk%u_matrix_opt(j,j2,kkk)
                 endif
              enddo
           else
              if (inwp%WFrot) then
                 do j2=1, chk%num_wann
                    local_contribution = local_contribution+ &
                         rotmat(j2,inwp%WFidx)*chk%u_matrix(j,j2,kkk)
                 enddo
              else
                 local_contribution = chk%u_matrix(j,inwp%WFidx,kkk)
              endif
           endif
           uu(j)=local_contribution
        enddo

        if(IERR.gt.0) KKKABS = KKK
        Npw = Nmatk - NLO

        ! Embiggen Nmatk-dependent variables if necessary.
        if (Nmatk > Nmat) then
           Nmat = Nmatk

           if (allocated(kx)) deallocate(kx,ky,kz, fj,dfj, BK,coef,cf)

           allocate ( kx  (Nmat), ky(Nmat), kz(  Nmat))
           allocate ( coef(Nmat), cf(Nmat), BK(3,Nmat) )
           allocate ( FJ(0:Lmax7, Nmat, NMT), DFJ(0:Lmax7, Nmat, NMT) )
        end if

        ! read in PW and local orbital wave vectors K
        read(unit_vector) (kx(ik), ky(ik), kz(ik), ik=1,Nmatk)
        !
        write(unit_out, "(' k-point',I4,' : ',A10,&
             &' = (',F8.4,',',F8.4,',',F8.4,')'/' ',55('-'))")&
             kkkabs, bname, fk

        ! transform K+k into global Cartesian coordinates
        do ik=1,Nmatk
           BK(:, ik) = matmul(stru%conv_rec, &
                &             (/ kx(ik)+fk(1), ky(ik)+fk(2), kz(ik)+fk(3) /))
        end do

        ! transform k into primitive fractional coordinates
        if(.not. stru%prim)then
           R  = matmul(stru%conv_rec, fk)
           fk = matmul(transpose(stru%prim_dir), R)
        else
           fk = TAU * fk
        endif

! load spherical Bessel functions and their derivatives
! -----------------------------------------------------------------
! FJ (l,:,Rmt) :      j_l(|K+k|*r) at r = Rmt_a for all PW's K+k
! DFJ(l,:,Rmt) : d/dr j_l(|K+k|*r) at r = Rmt_a for all PW's K+k
! -----------------------------------------------------------------
        call bessel(Npw, BK, NMT, rad, FJ, DFJ)

! read in individual eigen states
        COEF=0
        do I=1,NE
           read(unit_vector,IOSTAT=IERR) IE,E,IEABS
           read(unit_vector) CF(1:Nmatk)
           if (ie >= inwf%bmin .and. ie <= inwf%bmax) then
              coef(1:Nmatk) = coef(1:Nmatk) + uu(ie-inwf%bmin+1)*cf(1:Nmatk)
           end if
        end do

        ! normalize the eigen functions properly
        coef = prefac*coef

! set up the eigen state's augmentation coefficients
! ALM, BLM, and CLM for each atom in the unit cell
! ---------------------------------------------------------------------------
! in the muffin thin sphere around R + R_a one has:
!
! psi(r) = e^ikR Sum(lm) w_lm,a(|r-R-R_a|) Y*(*)_lm(T_a^-1(r-R-R_a))   with
!
! w_lm,a(r) = 4pi*i^l [ A_l,m,a *      u_l,a(r,E _l,a) +
!                       B_l,m,a * d/dE u_l,a(r,E _l,a) +
!                       C_l,m,a *      u_l,a(r,E'_l,a) ] * Rmt_a^2
!
! Here (*) stands for an optional additional complex conjugation on Y*_lm(…)
! WIEN95 : (*) =     and hence   *(*) = *
! WIEN97 : (*) = *   and hence   *(*) =
! ---------------------------------------------------------------------------
        il1 = npw
        latom  = 0
        aug_neq: do jneq=1,stru%Nneq
           imax = stru%Npt(jneq)
           aug_mult: do jmult=1,stru%mult(jneq)
              latom = latom + 1
              call augpw(stru, latom, rotloc_sym, pos_sym, Npw, &
                   &     Alm, Y, bk, coef, Nmatk)

              il=il1
              call auglo(stru, latom, rotloc_sym, pos_sym, il, &
                   &     Alm, Y, bk, coef, Nmatk)

              if (jmult == stru%mult(jneq)) il1=il

! add 4pi*i^l Rmt^2 factor to ALM, BLM and CLM and
! set up the eigen state's augmentation functions
! ---------------------------------------------------------------------------
! w_lm,a(r) = 4pi*i^l Rmt_a^2 * [ A_l,m,a *      u_l,a(r,E _l,a) +
!                                 B_l,m,a * d/dE u_l,a(r,E _l,a) +
!                                 C_l,m,a *      u_l,a(r,E'_l,a) ]
! ---------------------------------------------------------------------------
              PHS = 2*TAU * stru%RMT(jneq)**2
              LM = 0
              aug_l: do L=0,LMAX7
                 aug_m: do M=-L,L
                    LM=LM+1
                    do I=1,IMAX
                       ALMSUM=0
                       do irf=1,nrf
                          ALMSUM=ALMSUM+ALM(LM,irf)*RRAD(I,L,irf,JNEQ)
                       enddo
                       AUG(I,LM,LATOM)=ALMSUM*PHS
                    end do
                 end do aug_m
                 PHS = (0, 1._DPk)*PHS
              end do aug_l
           end do aug_mult
        end do aug_neq

        ! now evaluate the wave function on the grid
        do ig=1,Npt
           if(ireg(ig) == 0) then
              ! grid point in interstitial
              if(inwp%large)then
                 call wavint(inwp%rpts(1,ig),Npw,psi(ig),bk,coef,nmatk)
              else
                 psi(ig ) = 0
              endif
           else
              ! grid point in atomic sphere at R(ILAT) + R_a(IREG)
              ! load calc. Bloch factor e^ikR
              arg = dot_product(fk, ilat(:, ig))
              phs = cmplx(cos(arg), sin(arg), DPk)
              call wavsph(inwp%rpts(1,ig), phs, ireg(ig), iri(ig), &
                   &      psi(ig), Y, stru)
           endif
        end do

        tot = tot + psi
     end do kpoints
  end do vecfiles

  ! DEPHAS: A perfectly cromulent phase cannot always be achieved, but
  ! we will do our best.
  if (inwp%dephas) then
     PHSAV = 0
     do IG=1,Npt
        if(abs(tot(IG)) > PHASE_TOL) &
             PHSAV =  &
             PHSAV + mod(atan2(aimag(tot(ig)), real(tot(ig)))+TAU/2, TAU/2)
     end do
     PHSAV = PHSAV/real(Npt, DPk)
     PHS   = cmplx(cos(PHSAV), -sin(PHSAV), DPk)

     tot = tot * PHS
  end if

  write(unit_psink, 3020) abs(tot/kkk)**2
  write(unit_psiarg,3020) atan2(aimag(tot), real(tot))
3020 format(1P,10E16.8)
end subroutine main

subroutine wfopen(unit, file, iwf, iostat, status, form)
  use const, only: BUFSZ

  implicit none

  integer,      intent(in)            :: unit
  character(*), intent(in)            :: file
  integer,      intent(in)            :: iwf
  integer,      intent(out), optional :: iostat
  character(*), intent(in),  optional :: status
  character(*), intent(in),  optional :: form

  character(len=BUFSZ) :: tmpnam
  character(len=   11) :: f
  character(len=    7) :: s
  integer              :: idot

  f = 'FORMATTED'
  s = 'UNKNOWN'

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

pure logical function orthogonal(A)
  real(DPk), intent(in) :: A(:,:)

  orthogonal = all(abs( matmul(A, transpose(A)) - Id(size(A,1)) ) &
       &           < MAT_TOL)
end function orthogonal

pure logical function unitary(A)
  complex(DPk), intent(in) :: A(:,:)

  unitary = all(abs( matmul(A, conjg(transpose(A))) - Id(size(A,1))) &
       &        < MAT_TOL)
end function unitary

end module main_m


!!/---
!! Local Variables:
!! mode: f90
!! End:
!!\---
!!
!! Time-stamp: <2016-12-03 13:37:15 assman@faepop71.tu-graz.ac.at>
