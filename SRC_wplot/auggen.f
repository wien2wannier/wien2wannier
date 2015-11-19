!!! wien2wannier/SRC_wplot/auggen.f
!!!
!!! $Id: auggen.f 166 2014-02-03 09:39:24Z assmann $

      SUBROUTINE AUGGEN(REL,NAT,WHPSI)
      use struct
      use radgrd
      use lolog
      use loabc
      use atspdt
      use radfu
      use param
      use const, only: R8
!
      IMPLICIT REAL(R8) (A-H,O-Z)
!
      CHARACTER(5) WHPSI
      LOGICAL   REL, LARGE, rlo(1:nloat,0:lomax)
!
      DIMENSION AE(NRAD),BE(NRAD),VR(NRAD), &
                E(0:LMAX7),ELO(0:LOMAX,NLOAT),PEI(0:LMAX7)
      DIMENSION RAD1(NRAD,0:LMAX7,NRF),RAD2(NRAD,0:LMAX7,NRF)
      COMMON /WORK1/   A(NRAD),B(NRAD)
!---------------------------------------------------------------------
      CFAC = 1.0D0 / 274.074D0
      IF(.NOT.REL) CFAC = 1.0D-11
      LARGE = WHPSI .EQ. 'LARGE'
      IF(LARGE) CFAC = 1.0D0
!
!     << skip header of *.vsp file >>
      READ(unit_vsp,1000)

      WRITE(unit_out,2000)
      NLO = 0
      ALO = 0
      ILO = 0
      DO 10 JATOM=1,NAT
        IMAX = JRI(JATOM)
        ZNUC1 = Znuc (JATOM)
        RMT2 = RMT(JATOM)*RMT(JATOM)
        WRITE(unit_out,2010)JATOM,RMT(JATOM)
!
!       << read atomic potential r * V(r) and convert into Rydberg >>
        READ(unit_vsp,1010)
        READ(unit_vsp,1020)(VR(J),J=1,IMAX)
        READ(unit_vsp,1030)
        DO 20 J=1,IMAX
          VR(J) = 0.5D0 * VR(J)
   20   CONTINUE
!
!       << read augmentation energies and check for local orbitals >>
        READ(unit_vector) E
        READ(unit_vector) ELO
        DO L=0,LMAX7
         LAPW(L,JATOM)=.TRUE.
         IF (E(L).gt.150.) THEN
          E(L)=E(L)-200.d0
          LAPW(L,JATOM)=.FALSE.
         ENDIF
        ENDDO

        DO 30 L=0,LOMAX
         DO k=1,nloat
          rlo(k,l)=.false.
          IF(ELO(L,K).LT.995D0)THEN
            ilo(L,jatom)=ilo(L,jatom)+1
            IF(.NOT.lapw(l,jatom).AND.k.EQ.1) GOTO 666
            IF(k.EQ.nloat) THEN
             rlo(ilo(l,jatom),l)=.TRUE.
             GOTO 666
            ENDIF
 666        CONTINUE
            NLO=NLO+(2*L+1)*MULT(JATOM)
          ENDIF
         ENDDO
   30   CONTINUE
!
!       << regular augmentation functions >>
!       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        WRITE(unit_out,2020)
        DELE=2.0D-3
        DELEI=0.25D0/DELE
        DO 40 l=0,LMAX7
          FL=L
          EI=E(L)/2.0d0
!         << compute d/dE u_l(r,E) by finite differences >>
          E1=EI-DELE
          CALL OUTWIN(REL,VR,RM(1,JATOM),DX(JATOM),IMAX,E1, &
                      FL,UVB,DUVB,NODES,ZNUC1)
          CALL RINT13(REL,A,B,A,B,RNORM,JATOM)
          RNORM = 1.0D0/SQRT(RNORM)
          DO 50 M=1,IMAX
            AE(M) = RNORM * A(M)
            BE(M) = RNORM * B(M)
   50     CONTINUE
          UVB  = RNORM * UVB
          DUVB = RNORM * DUVB
          E1=EI+DELE
          CALL OUTWIN(REL,VR,RM(1,JATOM),DX(JATOM),IMAX,E1, &
                      FL,UVE,DUVE,NODES,ZNUC1)
          CALL RINT13(REL,A,B,A,B,RNORM,JATOM)
          RNORM = 1.0D0/SQRT(RNORM)
          UVE  = DELEI*(RNORM*UVE -UVB )
          DUVE = DELEI*(RNORM*DUVE-DUVB)
          DO 60 M=1,IMAX
            AE(M) = DELEI*(RNORM*A(M)-AE(M))
            BE(M) = DELEI*(RNORM*B(M)-BE(M))
   60     CONTINUE
!
!         << now compute u_l(r,E) >>
!
          CALL OUTWIN(REL,VR,RM(1,JATOM),DX(JATOM),IMAX,EI, &
                      FL,UV,DUV,NODES,ZNUC1)
          CALL RINT13(REL,A,B,A,B,RNORM,JATOM)
          RNORM = 1.0D0/SQRT(RNORM)
          DO 70 M=1,IMAX
            A(M) = RNORM*A(M)
            B(M) = RNORM*B(M)
   70     CONTINUE
          P(L,1,JATOM) = RNORM*UV
          DP(L,1,JATOM) = RNORM*DUV
!
!         << insure orthogonality of d/dE u_l(r,E) on u_l(r,E) >>
!
          CALL RINT13(REL,A,B,AE,BE,CROSS,JATOM)
          DO 80 M=1,IMAX
            AE(M) = (AE(M)-CROSS*A(M))
            BE(M) = (BE(M)-CROSS*B(M))
   80     CONTINUE
          P(L,2,JATOM) = UVE -CROSS*P(L,1,JATOM)
          DP(L,2,JATOM) = DUVE-CROSS*DP(L,1,JATOM)
          DO I=1,IMAX
            RAD1(I,L,1) = A(I)
            RAD1(I,L,2) = AE(I)
            RAD2(I,L,1) = B(I)
            RAD2(I,L,2) = BE(I)
          ENDDO
          CALL RINT13(REL,AE,BE,AE,BE,PEI(L),JATOM)
          WRITE(unit_out,2030) L,E(L),P(L,1,JATOM),DP(L,1,JATOM),P(L,2,JATOM),DP(L,2,JATOM)
   40   CONTINUE
!
!       << local orbitals >>
        DO 100 L=0,LOMAX
         irf=2
         DO 110 jlo=1,ilo(L,jatom)
          IF (lapw(l,jatom).or.(jlo.gt.1)) THEN
           irf=irf+1
           DELE=2.0D-3
           DELEI=0.25D0/DELE
           FL=L
           EI=elo(l,jlo)/2.d0
           IF(rlo(jlo,l)) THEN
            ei=elo(l,nloat)/2.d0
            kappa=l
            CALL diracout(rel,vr,RM(1,JATOM),dx(jatom),IMAX,    &
                          ei,fl,kappa,uv,duv,nodes,znuc1)
            CALL dergl(a,b,RM(1,jatom),dx(jatom),IMAX)
            DO m=1,IMAX
             r_m=RM(1,JATOM)*exp(dx(jatom)*(m-1))
             b(m)=b(m)*r_m/(2.d0*clight+(elo(l,jlo)- &
                   2.d0*vr(m)/r_m)/(2.d0*clight))
             b(m)=b(m)*clight
            ENDDO
           ELSE
            CALL outwin(rel,vr,RM(1,JATOM),dx(jatom),IMAX,   &
                        ei,fl,uv,duv,nodes,znuc1)
           ENDIF
           CALL RINT13(REL,A,B,A,B,RNORM,JATOM)
           RNORM = 1.0d0/SQRT(RNORM)
           DO M=1,IMAX
            RAD1(M,L,irf) = RNORM*A(M)
            RAD2(M,L,irf) = RNORM*B(M)
           ENDDO
           P(L,irf,jatom)  = RNORM*UV
           DP(L,irf,jatom) = RNORM*DUV
           CALL RINT13(REL,RAD1(1,L,1),RAD2(1,L,1),RAD1(1,L,irf),RAD2(1,L,irf),PI12LO,JATOM)
           CALL RINT13(REL,RAD1(1,L,2),RAD2(1,L,2),RAD1(1,L,irf),RAD2(1,L,irf),PE12LO,JATOM)
          ENDIF
        
          IF (LAPW(L,JATOM)) THEN
           XAC=(P(L,irf,jatom)*DP(L,2,JATOM)-DP(L,irf,jatom)*P (L,2,JATOM))*RMT2
           XBC=(P(L,1,jatom)*DP(L,irf,JATOM)-P(L,irf,JATOM)*DP(L,1,jatom))*RMT2
           XCC=XAC*(XAC+2.0D0*PI12LO) &
           + XBC*(XBC*PEI(L)+2.0D0*PE12LO)+1.0D0
           ALO(L,jlo,irf,JATOM) = 1.0D0/MAX( SQRT(XCC) , 0.005D0 )
           ALO(L,jlo,1,JATOM)=XAC*ALO(L,jlo,irf,JATOM)
           ALO(L,jlo,2,JATOM)=XBC*ALO(L,jlo,irf,JATOM)
          ELSE
           if (jlo.eq.1) then
            alonorm=sqrt(1.d0+(P(L,1,JATOM)/P(L,2,JATOM))**2*PEI(L))
            ALO(L,jlo,1,JATOM)=1.d0/alonorm
            ALO(L,jlo,2,JATOM)=-P(L,1,JATOM)/P(L,2,JATOM)/alonorm
           else
            xbc=-P(l,1,jatom)/P(L,irf,jatom)
            xac=sqrt(1+xbc**2+2*xbc*PI12LO)
            ALO(l,jlo,1,jatom)=1.d0/xac
            ALO(l,jlo,irf,jatom)=xbc/xac
           endif  
          ENDIF
  110    CONTINUE
  100   CONTINUE
!
        DO L=0,LOMAX
         do JLO=1,ilo(l,jatom)
          WRITE(unit_out,2070)L,(alo(l,jlo,jrf,jatom),jrf=1,nrf)
         enddo
        ENDDO

!   
!       << re-scale radial augmentation functions  >>
!       << with VR just being a working array here >>
        DO 135 I=1,IMAX
          VR(I) = CFAC / RM(I,JATOM)
  135   CONTINUE
        DO L=0,LMAX7
         IF(LARGE)THEN
          DO irf=1,nrf
           DO I=1,IMAX
            RRAD(I,L,irf,JATOM) = RAD1(I,L,irf) * VR(I)
           ENDDO
          ENDDO
         ELSE
          DO irf=1,nrf
           DO I=1,IMAX
            RRAD(I,L,irf,JATOM) = RAD2(I,L,irf) * VR(I)
           ENDDO
          ENDDO
         ENDIF
        ENDDO

   10 CONTINUE
      RETURN
!
 1000 FORMAT(//)
 1010 FORMAT(/////)
 1020 FORMAT(3X,4E19.12)
 1030 FORMAT(/////)
!
 2000 FORMAT(/' AUGMENTATION FUNCTIONS' &
             /' ----------------------' &
             /' psi_lm(r) = u_l(|r|) * Y_lm(r/|r|)')
 2010 FORMAT(/' augmentation functions for atom',i4,' (Rmt=',F6.3,')')
 2020 FORMAT(/' regular augmentation functions at r = Rmt' &
             /'  L    E(L)      u(r,E)     u''(r,E)   dE u(r,E)' &
             ,'  dE u''(r,E)  <u|dE u>  |dE u|^2')
 2030 FORMAT(I3,F8.3,1P,4E12.4,2E10.2)
 2070 FORMAT ('LO COEFFICIENT: l,A,B,C  ',i2,5X,6F12.5)
!
      END

!!/---
!! Local Variables:
!! mode: f90
!! End:
!!\---
