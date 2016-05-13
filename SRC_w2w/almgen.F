!!! wien2wannier/SRC_w2w/almgen.F

     SUBROUTINE almgen(alm,jatom,lfirst,NB,kkk)
     USE param
     USE struct
     USE xa3
     use xa
     use lolog, only: nlo,nlov,nlon,ilo,lapw,n_rad
     use gener, only: br1, br2

     IMPLICIT REAL(R8) (A-H,O-Z)

     real(r8)     h_al(iblock),h_bl(iblock)

     COMPLEX(C16)       YL((LMAX2+1)*(LMAX2+1))
     COMPLEX(C16)       PHSHEL,CFAC
     complex(c16) h_yl((LMAX2+1)*(LMAX2+1),iblock)
     complex(c16) h_alyl((LMAX2+1)*(LMAX2+1),iblock)
     complex(c16) h_blyl((LMAX2+1)*(LMAX2+1),iblock)
     complex(c16) ALMt((LMAX2+1)*(LMAX2+1),NB,NRF,MULT(JATOM))
     complex(c16) ALM(NB,NRF,(LMAX2+1)*(LMAX2+1),MULT(JATOM))

     COMMON /ATSPDT/  P(0:LMAX2,nrf),DP(0:LMAX2,nrf)
     COMMON /loabc/   alo(0:lomax,nloat,nrf)

     TWOPI=2.D0*PI
     FAC=4.0D0*PI*RMT(JATOM)**2/SQRT(VOL)

     N=size(kkk)
     ALMt=0
     DO I=1,N
      BKX(I)=XK(kkk)+GX(I,kkk)                    ! plane wave Q+k
      BKY(I)=YK(kkk)+GY(I,kkk)
      BKZ(I)=ZK(kkk)+GZ(I,kkk)
     ENDDO
     CALL HARMON(N,BKX,BKY,BKZ,LMAX2,FJ,DFJ,RMT(JATOM))    ! calculate Bessel functions j_l(|Q+k|*RMT)
     DO 777 MU=1,MULT(JATOM)    !  loop over equivalent atoms
      latom=lfirst-1+mu
      DO 120 ii=1,N-(nlo+nlon+nlov),iblock                     ! matching of plane wave and atomic functions
       i3=0
       DO 121 i=ii,min(ii+iblock-1,N-(nlo+nlon+nlov))
        i3=i3+1
        BK(1)=BKX(I)*BR1(1,1)+BKY(I)*BR1(1,2)+BKZ(I)*BR1(1,3)   ! transform Q+k to cartesian coordinates
        BK(2)=BKX(I)*BR1(2,1)+BKY(I)*BR1(2,2)+BKZ(I)*BR1(2,3)
        BK(3)=BKX(I)*BR1(3,1)+BKY(I)*BR1(3,2)+BKZ(I)*BR1(3,3)
        CALL YLM (BK,LMAX2,YL)                                  ! calculate Y_lm(Q+k)
        ARG1=BKX(I)*POS(1,LATOM)*TWOPI
        ARG2=BKY(I)*POS(2,LATOM)*TWOPI
        ARG3=BKZ(I)*POS(3,LATOM)*TWOPI
        PHSHEL=EXP((0,1)*(ARG1+ARG2+ARG3))                       ! exp(i(Q+k)R_j)
        INDEX=0
        DO  L=0,LMAX2
         MAX=2*L+1
         DO  M=1,MAX
          INDEX=INDEX+1
          h_yl(index,i3)=conjg(yl(index))*phshel
         ENDDO
        ENDDO
  121  CONTINUE
       INDEX=0
       DO  L=0,LMAX2
        i3=0
        do  i=ii,min(ii+iblock-1,N-(nlo+nlon+nlov))
         i3=i3+1
         if (lapw(l)) then
          h_AL(i3)=DFJ(L,I)*P(L,2)-FJ(L,I)*DP(L,2)
          h_BL(i3)=FJ(L,I)*DP(L,1)-DFJ(L,I)*P(L,1)
         else
          h_AL(i3)=FJ(L,I)/P(L,1)/RMT(JATOM)**2
          h_BL(i3) = 0.d0
         end if
        enddo
        MAX=2*L+1
        DO M=1,MAX
         INDEX=INDEX+1
         i3=0
         do  i=ii,min(ii+iblock-1,N-(nlo+nlon+nlov))
          i3=i3+1
          h_alyl(index,i3)=h_AL(i3)*h_YL(INDEX,i3)
          h_blyl(index,i3)=h_BL(i3)*h_YL(INDEX,i3)
         enddo
        ENDDO
       ENDDO
       ibb=min(iblock,N-(nlo+nlon+nlov)-ii+1)
#ifndef _COMPLEX_
       lda=2*(LMAX2+1)*(LMAX2+1)
#else
       lda=(LMAX2+1)*(LMAX2+1)
#endif
           ldc=lda
           ldb=nmat
#ifndef _COMPLEX_
           call dgemm('N','N',2*index,nb,ibb,1.d0, &
                h_alyl,lda,a(ii,1,kkk),ldb,1.d0, &
                almt(1,1,1,mu),ldc)
           call dgemm('N','N',2*index,nb,ibb,1.d0, &
                h_blyl,lda,a(ii,1,kkk),ldb,1.d0, &
                almt(1,1,2,mu),ldc)
#else
           call zgemm('N','N',index,nb,ibb,(1.d0,0.d0), &
                h_alyl,lda,a(ii,1,kkk),ldb,(1.d0,0.d0), &
                almt(1,1,1,mu),ldc)
           call zgemm('N','N',index,nb,ibb,(1.d0,0.d0), &
                h_blyl,lda,a(ii,1,kkk),ldb,(1.d0,0.d0), &
                almt(1,1,2,mu),ldc)
#endif
  120 CONTINUE
      IF (nlo.ne.0) THEN
       i=n-(nlo+nlon)
       DO L=0,LoMAX
        mirf=n_rad(l)
        DO jlo=1,ilo(l)
         DO jneq=1,mult(jatom)
          DO M1=-l,l
           i=i+1
           BK(1)=BKX(I)*BR1(1,1)+BKY(I)*BR1(1,2)+BKZ(I)*BR1(1,3)
           BK(2)=BKX(I)*BR1(2,1)+BKY(I)*BR1(2,2)+BKZ(I)*BR1(2,3)
           BK(3)=BKX(I)*BR1(3,1)+BKY(I)*BR1(3,2)+BKZ(I)*BR1(3,3)
           CALL YLM (BK,LOMAX,YL)
           ARG=BKX(I)*POS(1,LATOM)+BKY(I)*POS(2,LATOM)+BKZ(I)*POS(3,LATOM)
           PHSHEL=EXP((0,1)*(ARG*TWOPI))
           DO NUM=1,NB
            PHS(NUM)=PHSHEL*A(I,NUM,kkk)
           ENDDO
           DO M=-l,l
            index=l*(l+1)+m+1
            DO NUM=1,NB
             DO irf=1,mirf
              ALMt(index,num,irf,mu)=ALMt(INDEX,num,irf,mu)+ &
              ALo(l,jlo,irf)*dconjg(YL(INDEX))*PHS(NUM)
             ENDDO
            ENDDO
           ENDDO
          ENDDO
         ENDDO
        ENDDO
       ENDDO
      ENDIF

      INDEX=0
      DO L=0,LMAX2
       mirf=n_rad(l)
       MAX=2*L+1
       CFAC=(0,1)**L
       DO M=1,MAX
        INDEX=INDEX+1
        DO irf=1,mirf
         DO NUM=1,NB
          ALM(num,irf,index,mu)=ALMt(INDEX,NUM,irf,mu)*FAC*CFAC
         ENDDO
        ENDDO
       ENDDO
      ENDDO
 777 CONTINUE
     RETURN
     END


!!/---
!! Local Variables:
!! mode: f90
!! End:
!!\---
!!
!! Time-stamp: <2016-04-08 17:11:36 assman@faepop36.tu-graz.ac.at>
