!!! wien2wannier/SRC_wplot/auglo.f
!!!
!!! $Id: auglo.f 166 2014-02-03 09:39:24Z assmann $

      SUBROUTINE AUGLO(LATOM,IL,ALM,ROTLOC,Y,bk,coef,nmat)
      use struct
      use lolog
      use loabc
      use param
      use const, only: R8, C16
!
      IMPLICIT REAL(R8) (A-H,O-Z)
      COMPLEX(C16) ALM((LMAX7+1)*(LMAX7+1),NRF)
      
      COMPLEX(C16) COEF(nmat) !changed by pwissgott
      real(r8) BK(3,NMAT)
      COMPLEX(C16) PHS,PHSSUM(-LOMAX:LOMAX),Y((LOMAX+1)*(LOMAX+1))
      DIMENSION  RK(3),ROTLOC(3,3,*)
!
      JATOM = IATNR(LATOM)
      DO 10 L=0,LOMAX
         DO 20 jlo=1,ilo(L,JATOM)
            PHSSUM = (0.0D0,0.0D0)
            DO 30 JNEQ=1,MULT(JATOM)
               DO M1=-L,L
                  IL=IL+1
                  DO I=1,3
                     RK(I)=ROTLOC(I,1,LATOM)*BK(1,IL)+  &
                     ROTLOC(I,2,LATOM)*BK(2,IL)+ &
                     ROTLOC(I,3,LATOM)*BK(3,IL)
                  ENDDO
                  CALL YLM(RK,LOMAX,Y)
                  ARG=BK(1,IL)*POS(1,LATOM)+BK(2,IL)*POS(2,LATOM) &
                  + BK(3,IL)*POS(3,LATOM)
                  PHS=COEF(IL)*DCMPLX(COS(ARG),SIN(ARG))
                  LM = L*L
                  DO M=-L,L
                     LM = LM + 1
                     PHSSUM(M) = PHSSUM(M) + PHS * DCONJG( Y(LM) )
                  ENDDO
               ENDDO
 30         CONTINUE
            LM = L*L
            DO M=-L,L
               LM=LM+1
               DO irf=1,NRF     !changed be p.wissgott
                  ALM(LM,irf)=ALM(LM,irf)+PHSSUM(M)*ALO(L,jlo,irf,JATOM)
               ENDDO
            ENDDO
 20      CONTINUE
 10   CONTINUE
      RETURN
      END

!!/---
!! Local Variables:
!! mode: fortran
!! End:
!!\---
