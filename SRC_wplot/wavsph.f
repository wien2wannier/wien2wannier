!!! wien2wannier/SRC_wplot/wavsph.f

      SUBROUTINE WAVSPH(R,BFAC,IAT,IR,PSI,Y)
      use struct
      use radgrd                     
      use work                     
      use param

      IMPLICIT REAL(R8) (A-H,O-Z)
      DIMENSION  R(3)
      COMPLEX(C16) BFAC,PSI
!
! evaluation of the wave function in the a muffin tin sphere around R+R_a
!
! psi(r) = e^ikR Sum(lm) w_lm,a(|r-R-R_a|) Y*(*)_lm(T_a^-1(r-R-R_a))
!
! Here (*) stands for an optional additional complex conjugation on Y*_lm(...)
! WIEN95 : (*) =     and hence   *(*) = *
! WIEN97 : (*) = *   and hence   *(*) = 
! -----------------------------------------------------------------------
! Input:
! R    -- grid point in the local coordinate system of atom a
!         i.e the value s_a := T_a^-1(r-R-R_a)
! BFAC -- the Bloch factor e^ikR
! IAT  -- the atom a within the unit cell
! IR   -- the radial interval [r_i,r_i+1] the value |r-R-R_a| falls in
!
! COMMON /WORK  /
! AUG(:,lm,a) -- the augmentation functions w_lm,a(r) on the radial mesh
!
! COMMON /STRUCT/  structural information
! COMMON /RADGRD/  radial grid information
!
! Output:
! PSI  -- the wave function psi(r)
!
! Working arrays:
! Y    -- to hold Y_lm(...)
! -----------------------------------------------------------------------
!
      COMPLEX(C16) PHS,PHSLM,Y((LMAX7+1)*(LMAX7+1))
!
!     << Y*(*)_lm(T_a^-1(r-R-R_a)) for all lm >>
      CALL YLM(R,LMAX7,Y)
!:17[
      IF(.NOT.KCONJG)THEN
!       << WIEN95 convention : (*) =   >>
        DO 5 LM=1,(LMAX7+1)*(LMAX7+1)
          Y(LM) = DCONJG(Y(LM))
   5    CONTINUE
      ENDIF
!:17]
!
!     << prepare radial interpolation >>
! --------------------------------------------------------------------------
! r in [r1,r2] : w(r) = [ w(r1) ( r - r2) + w(r2) ( r1 - r ) ] / ( r2 - r1 )
! --------------------------------------------------------------------------
      JATOM = IATNR(IAT)
      RR = SQRT( R(1)*R(1) + R(2)*R(2) + R(3)*R(3) )
      R1 = RM(IR,JATOM)
      R2 = RM(IR+1,JATOM)
      W1 = (R2-RR)/(R2-R1)
      W2 = (RR-R1)/(R2-R1)
!
!     << Sum(lm) ... >>
      PHS = (0.0D0,0.0D0)
      DO 10 LM=1,(LMAX7+1)*(LMAX7+1)
        PHSLM = W1*AUG(IR,LM,IAT) + W2*AUG(IR+1,LM,IAT)
        PHS = PHS + PHSLM * Y(LM)
   10 CONTINUE
!
!     << psi(r) >>
      PSI = BFAC * PHS
!
      RETURN
      END


!!/---
!! Local Variables:
!! mode: f90
!! End:
!!\---
!!
!! Time-stamp: <2015-05-23 19:58:48 elias>
