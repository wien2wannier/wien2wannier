!!! wien2wannier/SRC_wplot/augpw.f

      subroutine AUGPW(latom, Npw, Alm, rotloc, Y, bk, coef, Nmat, iatnr)
      use struct, only: POS, RMT
      use atspdt, only: P, DP
      use bessfu, only: irad, fj, dfj
      use lolog,  only: lapw
      use param,  only: DPk, kconjg, Lmax7, Nrf

      implicit none

      integer,      intent(in)  :: Latom, Npw, Nmat, iatnr(*)
      real(DPk),    intent(in)  :: rotloc(3,3,*), BK(3,Nmat)
      complex(DPk), intent(out) :: ALM((LMAX7+1)*(LMAX7+1),nrf)
      complex(DPk), intent(out) :: Y((LMAX7+1)*(LMAX7+1))
      complex(DPk), intent(in)  :: COEF(nmat)

      complex(DPk) :: PHS, PHSLM
      real(DPk)    :: RK(3), arg, al, bl
      integer      :: i, imt, jatom, iPW, lm, l, m

! The PW part of the augmentation coefficients A_lm,a and B_lm,a
! of a given eigen state at a given atom a in the unit cell
! ----------------------------------------------------------------------------
! Input:
! LATOM    -- the atom a
! NPW      -- current number of PW basis functions
! ROTLOC   -- the local rotation matrices T_a of each atom
!             [ stored as (T_a^-1)_ij = ROTLOC(i,j,a) ]
!
! Output:
! ALM  -- the PW part of the augmentation coefficients for each (l,m)
!
! Working arrays:
! Y        -- to hold Y_lm(...)
! ----------------------------------------------------------------------------
! X_l,m,a = Sum(K) c_K/sqrt(V) exp(i(K+k)R_a) Y(*)_lm(T_a^-1(K+k)) X_l,a(|K+k|)
!
! Here (*) stands for an optional complex conjugation on Y_lm(...)
! WIEN95 : (*) =
! WIEN97 : (*) = *
!
! R_a   : center of atom a
! T_a   : local rotation matrix of atom a
! X_l,a : PW augmentation coefficients for atom a
!
! Note: Letting  R_a = Q_a(R_0) + t_a  yields
!
!       exp(i(K+k)R_a) = exp(i(K+k)t_a + i[Q_a^-1(K+k)]R_0)
!
!       which is the expression used in LAPW1 and LAPW2
! ---------------------------------------------------------------------------

      JATOM = iatnr(LATOM)
      IMT   = IRAD (JATOM)
!
!     << initialize ALM and BLM >>
        ALM = 0
!
      DO 20 IPW=1,NPW
!
!       << Y*(*)_lm(T_a^-1(K+k)) >>
        DO 30 I=1,3
          RK(I) = ROTLOC(I,1,LATOM)*BK(1,IPW) &
                + ROTLOC(I,2,LATOM)*BK(2,IPW) &
                + ROTLOC(I,3,LATOM)*BK(3,IPW)
   30   CONTINUE
        CALL YLM(RK,LMAX7,Y)

!:17[
        IF(.NOT.KCONJG)THEN
!         << WIEN95 convention : (*) =   >>
          DO 35 LM=1,(LMAX7+1)*(LMAX7+1)
            Y(LM) = conjg(Y(LM))
   35     CONTINUE
        ENDIF
!:17]
!
!       << c_K/sqrt(V) * exp(i(K+k)R_a) >>
        ARG = BK(1,IPW)*POS(1,LATOM) + BK(2,IPW)*POS(2,LATOM) &
            + BK(3,IPW)*POS(3,LATOM)
        PHS = COEF(IPW) * cmplx(cos(arg), sin(arg), DPk)
!
!       << A_l,a and B_l,a (without Rmt^2 factor) >>
! -----------------------------------------------------------------
! A_a,l = - [ d/dE d/dr u_l,a(Rmt,E)      j_l(Rmt*|K+k|) -
!             d/dE      u_l,a(Rmt,E) d/dr j_l(Rmt*|K+k|) ]
!
! B_a,l = - [           u_l,a(Rmt,E) d/dr j_l(Rmt*|K+k|) -
!                  d/dr u_l,a(Rmt,E)      j_l(Rmt*|K+k|) ]
! -----------------------------------------------------------------
        LM = 0
        DO 40 L=0,LMAX7
         if (lapw(l,jatom)) then
          AL = DFJ(L,IPW,IMT)* P(L,2,JATOM) - &
                FJ(L,IPW,IMT)*DP(L,2,JATOM)
          BL =  FJ(L,IPW,IMT)*DP (L,1,JATOM) -  &
               DFJ(L,IPW,IMT)* P (L,1,JATOM)
          else
           AL = FJ(L,IPW,IMT)/P(L,1,JATOM)/RMT(JATOM)**2
           BL= 0.d0
          endif
          DO 50 M=-L,L
            LM = LM + 1
            PHSLM = PHS * conjg( Y(LM) )
            ALM(LM,1) = ALM(LM,1) + PHSLM * AL
            ALM(LM,2) = ALM(LM,2) + PHSLM * BL
   50     CONTINUE
   40   CONTINUE
   20 CONTINUE
!
      END


!!/---
!! Local Variables:
!! mode: f90
!! End:
!!\---
!!
!! Time-stamp: <2015-12-29 19:57:06 assman@faepop36.tu-graz.ac.at>
