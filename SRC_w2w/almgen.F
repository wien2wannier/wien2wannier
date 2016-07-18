!!! wien2wannier/SRC_w2w/almgen.F

module     almgen_m; contains
subroutine almgen(Alm,jatom,lfirst,NB,kkk)
  use param,  only: Nmat, lomax, NRF, iblock, Lmax2
  use struct, only: RMT, vol, mult, pos
  use xa3,    only: vecsz, XK,YK,ZK, BKX,BKY,BKZ, GX,GY,GZ, A
  use xa,     only: fj, dfj, bk, phs
  use lolog,  only: nlo,nlov,nlon,ilo,lapw,n_rad
  use gener,  only: br1
  use atspdt, only: P, DP
  use loabc,  only: alo
  use const,  only: R8, C16, TAU

  !! procedure includes
  use harmon_m
  use Ylm_m

  implicit none

  integer,      intent(in)  :: jatom, lfirst, NB, kkk
  complex(C16), intent(out) :: Alm(NB,NRF,(LMAX2+1)*(LMAX2+1),MULT(JATOM))

  real(R8) :: h_al(iblock), h_bl(iblock)

  complex(C16) :: YL((LMAX2+1)*(LMAX2+1))
  complex(C16) :: PHSHEL,CFAC
  complex(C16) :: h_yl((LMAX2+1)*(LMAX2+1),iblock)
  complex(C16) :: h_alyl((LMAX2+1)*(LMAX2+1),iblock)
  complex(C16) :: h_blyl((LMAX2+1)*(LMAX2+1),iblock)
  complex(C16) :: ALMt((LMAX2+1)*(LMAX2+1),NB,NRF,MULT(JATOM))

  real(R8) :: fac, arg, arg1, arg2, arg3
  integer  :: N, num, i,ii,i3, latom, mu, index, l, max, m,m1, irf,mirf, ibb, lda,ldb,ldc, jlo,jneq

  FAC=2*TAU*RMT(JATOM)**2/sqrt(VOL)

  N=vecsz(kkk)
  ALMt=0
  do I=1,N
     BKX(I)=XK(kkk)+GX(I,kkk)                    ! plane wave Q+k
     BKY(I)=YK(kkk)+GY(I,kkk)
     BKZ(I)=ZK(kkk)+GZ(I,kkk)
  end do
  call harmon(N,BKX,BKY,BKZ,LMAX2,FJ,DFJ,RMT(JATOM))    ! calculate Bessel functions j_l(|Q+k|*RMT)
  do MU=1,MULT(JATOM)    !  loop over equivalent atoms
     latom=lfirst-1+mu
     do ii=1,N-(nlo+nlon+nlov),iblock                     ! matching of plane wave and atomic functions
        i3=0
        do i=ii,min(ii+iblock-1,N-(nlo+nlon+nlov))
           i3=i3+1
           BK(1)=BKX(I)*BR1(1,1)+BKY(I)*BR1(1,2)+BKZ(I)*BR1(1,3)   ! transform Q+k to cartesian coordinates
           BK(2)=BKX(I)*BR1(2,1)+BKY(I)*BR1(2,2)+BKZ(I)*BR1(2,3)
           BK(3)=BKX(I)*BR1(3,1)+BKY(I)*BR1(3,2)+BKZ(I)*BR1(3,3)
           call YLM (BK,LMAX2,YL)                                  ! calculate Y_lm(Q+k)

           ARG1=BKX(I)*POS(1,LATOM)*TAU
           ARG2=BKY(I)*POS(2,LATOM)*TAU
           ARG3=BKZ(I)*POS(3,LATOM)*TAU

           PHSHEL=EXP((0,1)*(ARG1+ARG2+ARG3))                       ! exp(i(Q+k)R_j)
           INDEX=0
           do L=0,LMAX2
              MAX=2*L+1
              do M=1,MAX
                 INDEX=INDEX+1
                 h_yl(index,i3)=conjg(yl(index))*phshel
              end do
           end do
        end do

        INDEX=0
        do L=0,LMAX2
           i3=0
           do i=ii,min(ii+iblock-1,N-(nlo+nlon+nlov))
              i3=i3+1
              if (lapw(l)) then
                 h_AL(i3)=DFJ(L,I)*P(L,2)-FJ(L,I)*DP(L,2)
                 h_BL(i3)=FJ(L,I)*DP(L,1)-DFJ(L,I)*P(L,1)
              else
                 h_AL(i3)=FJ(L,I)/P(L,1)/RMT(JATOM)**2
                 h_BL(i3) = 0.d0
              end if
           end do

           MAX=2*L+1
           do M=1,MAX
              INDEX=INDEX+1
              i3=0
              do i=ii,min(ii+iblock-1,N-(nlo+nlon+nlov))
                 i3=i3+1
                 h_alyl(index,i3)=h_AL(i3)*h_YL(INDEX,i3)
                 h_blyl(index,i3)=h_BL(i3)*h_YL(INDEX,i3)
              end do
           end do
        end do

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
     end do

     if (nlo.ne.0) then
        i=n-(nlo+nlon)
        do L=0,LoMAX
           mirf=n_rad(l)
           do jlo=1,ilo(l)
              do jneq=1,mult(jatom)
                 do M1=-l,l
                    i=i+1
                    BK(1)=BKX(I)*BR1(1,1)+BKY(I)*BR1(1,2)+BKZ(I)*BR1(1,3)
                    BK(2)=BKX(I)*BR1(2,1)+BKY(I)*BR1(2,2)+BKZ(I)*BR1(2,3)
                    BK(3)=BKX(I)*BR1(3,1)+BKY(I)*BR1(3,2)+BKZ(I)*BR1(3,3)
                    CALL YLM (BK,LOMAX,YL)
                    ARG=BKX(I)*POS(1,LATOM)+BKY(I)*POS(2,LATOM)+BKZ(I)*POS(3,LATOM)
                    PHSHEL=EXP((0,1)*(ARG*TAU))
                    DO NUM=1,NB
                       PHS(NUM)=PHSHEL*A(I,NUM,kkk)
                    ENDDO
                    do M=-l,l
                       index=l*(l+1)+m+1
                       do NUM=1,NB
                          do irf=1,mirf
                             ALMt(index,num,irf,mu)=ALMt(INDEX,num,irf,mu)+ &
                                  ALo(l,jlo,irf)*dconjg(YL(INDEX))*PHS(NUM)
                          end do
                       end do
                    end do
                 end do
              end do
           end do
        end do
     end if

     INDEX=0
     do L=0,LMAX2
        mirf=n_rad(l)
        MAX=2*L+1
        CFAC=(0,1)**L
        do M=1,MAX
           INDEX=INDEX+1
           do irf=1,mirf
              do NUM=1,NB
                 ALM(num,irf,index,mu)=ALMt(INDEX,NUM,irf,mu)*FAC*CFAC
              end do
           end do
        end do
     end do
  end do
end subroutine almgen
end module     almgen_m

!!/---
!! Local Variables:
!! mode: f90
!! End:
!!\---
!!
!! Time-stamp: <2016-07-18 15:47:49 assman@faepop71.tu-graz.ac.at>
