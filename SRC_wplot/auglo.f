!!! wien2wannier/SRC_wplot/auglo.f

      subroutine auglo(latom,il,Alm,rotloc,Y,bk,coef,nmat,stru)
      use struct,    only: POS
      use lolog,     only: iLO
      use loabc,     only: Alo
      use param,     only: DPk, Lmax7, LOmax, Nrf
      use structmod, only: struct_t

      implicit none

      integer,        intent(in)    :: latom, Nmat
      integer,        intent(inout) :: il
      complex(DPk),   intent(inout) :: Alm((Lmax7+1)*(Lmax7+1),Nrf)
      complex(DPk),   intent(in)    :: coef(Nmat)
      real(DPk),      intent(in)    :: BK(3,NMAT), rotloc(3,3,*)
      type(struct_t), intent(in)    :: stru

      COMPLEX(DPk) :: PHS, PHSSUM(-LOMAX:LOMAX), Y((LOMAX+1)*(LOMAX+1))
      real(DPk)    :: RK(3), arg
      integer      :: i, l, lm, m, m1, irf, jlo, jatom, jneq

      JATOM = stru%neq2at(LATOM)
      DO 10 L=0,LOMAX
         DO 20 jlo=1,ilo(L,JATOM)
            PHSSUM = (0.0D0,0.0D0)
            DO 30 JNEQ=1,stru%MULT(JATOM)
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
!! mode: f90
!! End:
!!\---
!!
!! Time-stamp: <2015-12-29 18:35:48 assman@faepop36.tu-graz.ac.at>
