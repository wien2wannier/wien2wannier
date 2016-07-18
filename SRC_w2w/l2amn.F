!!! wien2wannier/SRC_w2w/l2amn.F

module     l2amn_m; contains
subroutine l2amn(NB,NPROJ,NUM_KPTS)
  use param,   only: Lmax2, iblock, Ndif, NRF, Nmat, lomax, &
       unit_vsp, unit_out, unit_amn
  use const,   only: R8, C16, TAU
  use struct,  only: Nat, mult, rotloc, RMT, vol, rotij, transij, pos
  use xa,      only: BK, FJ, DFJ, PHS, BKrot, BKRloc
  use xa3,     only: A, XK,YK,ZK, GX,GY,GZ, BKX,BKY,BKZ, vecsz
  use amn_mmn, only: C
  use lolog,   only: nlo,nlov,nlon,loor,ilo,lapw
  use gener,   only: br1
  use loabc,   only: alo
  use atspdt,  only: P, DP
  use const,   only: R8, C16
  use radfu,   only: RF1, RF2

  !! procedure includes
  use rint13_m
  use rotate_m
  use harmon_m
  use atpar_m
  use Ylm_m

  implicit none

  integer, intent(in) :: NB, Nproj, num_kpts

  COMPLEX(C16) :: YL((LMAX2+1)*(LMAX2+1))        ! list of spherical harmics
  COMPLEX(C16) :: projection(NB,NPROJ,NUM_KPTS)
  COMPLEX(C16) :: PHSHEL,CFAC

  real(R8)     :: h_al(iblock),h_bl(iblock)
  complex(C16) :: h_yl((LMAX2+1)*(LMAX2+1),iblock)
  complex(C16) :: h_alyl((LMAX2+1)*(LMAX2+1),iblock)
  complex(C16) :: h_blyl((LMAX2+1)*(LMAX2+1),iblock)
  complex(C16) :: alm((LMAX2+1)*(LMAX2+1),nb,ndif,nrf)

  integer  :: i,ii,i3,j,k,l, m,mu,m1, num,num1,num2, n,n1,n2
  integer  :: jatom,jlo,jneq, index, ibb,iscf,irf, lda,ldb,ldc, latom,lfirst
  integer  :: kkk, max, itap,itape,jtape
  real(R8) :: fac, arg, pi12lo

  do k=1,NUM_KPTS
     do n2=1,NPROJ
        do n1=1,NB
           projection(n1,n2,k)=0
        enddo
     enddo
  enddo
!     ---------------------------------
!     START LOOP FOR ALL ATOMS BY DO 50
!     ---------------------------------
  rewind(unit_vsp)
  READ(unit_vsp, '(50X,I2,//)') ISCF
  LFIRST=1
  atoms: DO JATOM=1,NAT
     IF(JATOM.GT.1) LFIRST=LFIRST + MULT(JATOM-1)
     ITAP=30+JATOM
     itape=9
     jtape=unit_vsp
     rewind(itape)
     call atpar(JATOM, itape, jtape)
     rewind(itape)
     write(unit_out,*)'******************************************'
     do i=1,3
        write(unit_out, '(3f12.6)')(rotloc(i,j,jatom),j=1,3)
     enddo
     write(unit_out,*)

!................................
     kpoints: DO KKK=1,NUM_KPTS      !   loop over k
        N=vecsz(kkk)
        DO I=1,N
           BKX(I)=XK(kkk)+GX(I,kkk)                    ! plane wave Q+k
           BKY(I)=YK(kkk)+GY(I,kkk)
           BKZ(I)=ZK(kkk)+GZ(I,kkk)
        ENDDO
        call harmon(N,BKX,BKY,BKZ,LMAX2,FJ,DFJ,RMT(JATOM))    ! calculate Bessel functions j_l(|Q+k|*RMT)

        FAC=2*TAU*RMT(JATOM)**2/sqrt(VOL)
        ALM=0

        equiv1: DO MU=1,MULT(JATOM)    !  loop over equivalent atoms
           !  write(*,*)"debugatom",jatom,mu
           latom=lfirst-1+mu
           if (kkk.eq.1) then
              do i=1,3
                 write(unit_out, '(3f12.6,4x,f12.4)') &
                      (rotij(i,j,latom),j=1,3),transij(i,latom)
              enddo
              write(unit_out,*)
           endif

           DO ii=1,N-(nlo+nlon+nlov),iblock ! matching of plane wave and atomic functions
              i3=0
              DO i=ii,min(ii+iblock-1,N-(nlo+nlon+nlov))
                 i3=i3+1
                 bk(1)=bkx(i)
                 bk(2)=bky(i)
                 bk(3)=bkz(i)
                 call ROTATE (bk,rotij(1,1,latom),bkrot) ! go to coord. of latom
                 bk(1)=bkrot(1)*br1(1,1)+bkrot(2)*br1(1,2)+bkrot(3)*br1(1,3) ! transform Q+k to cartesian coordinates
                 bk(2)=bkrot(1)*br1(2,1)+bkrot(2)*br1(2,2)+bkrot(3)*br1(2,3)
                 bk(3)=bkrot(1)*br1(3,1)+bkrot(2)*br1(3,2)+bkrot(3)*br1(3,3)
                 call ROTATE (bk,rotloc(1,1,jatom),bkrloc) ! go to local (rotated) coords
                 CALL YLM (bkrloc,LMAX2,YL) ! calculate Y_lm(Q+k)
                 arg=bkrot(1)*pos(1,lfirst)+  &
                      bkrot(2)*pos(2,lfirst)+  &
                      bkrot(3)*pos(3,lfirst)+  &
                      bkx(i)*transij(1,latom)+bky(i)*transij(2,latom)+bkz(i)*transij(3,latom)
                 phshel=exp((0,1)*(arg*TAU))
                 INDEX=0
                 DO L=0,LMAX2
                    MAX=2*L+1
                    DO M=1,MAX
                       INDEX=INDEX+1
                       h_yl(index,i3)=conjg(yl(index))*phshel
                    ENDDO
                 ENDDO
              ENDDO
              INDEX=0
              DO L=0,LMAX2
                 i3=0
                 do i=ii,min(ii+iblock-1,N-(nlo+nlon+nlov))
                    i3=i3+1
                    if (lapw(l)) then
                       h_AL(i3)=DFJ(L,I)*P(L,2)-FJ(L,I)*DP(L,2)
                       h_BL(i3)=FJ(L,I)*DP(L,1)-DFJ(L,I)*P(L,1)
                    else
                       h_AL(i3)=FJ(L,I)/P(L,1)/RMT(JATOM)**2
                       h_BL(i3) = 0.d0
                    endif
                 enddo
                 MAX=2*L+1
                 DO M=1,MAX
                    INDEX=INDEX+1
                    i3=0
                    do i=ii,min(ii+iblock-1,N-(nlo+nlon+nlov))
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
                   alm(1,1,mu,1),ldc)
              call dgemm('N','N',2*index,nb,ibb,1.d0, &
                   h_blyl,lda,a(ii,1,kkk),ldb,1.d0, &
                   alm(1,1,mu,2),ldc)
#else
              call zgemm('N','N',index,nb,ibb,(1.d0,0.d0), &
                   h_alyl,lda,a(ii,1,kkk),ldb,(1.d0,0.d0), &
                   alm(1,1,mu,1),ldc)
              call zgemm('N','N',index,nb,ibb,(1.d0,0.d0), &
                   h_blyl,lda,a(ii,1,kkk),ldb,(1.d0,0.d0), &
                   alm(1,1,mu,2),ldc)
#endif
           ENDDO
!********  LOCAL ORBITALS *****************************
           lo: IF (nlo.ne.0) THEN
              i=n-(nlo+nlon)
              DO L=0,LoMAX
                 DO jlo=1,ilo(l)
                    DO jneq=1,mult(jatom)
                       DO M1=-l,l
                          i=i+1
                          BK(1)=BKX(I)
                          BK(2)=BKY(I)
                          BK(3)=BKZ(I)
                          CALL ROTATE (BK,ROTIJ(1,1,LATOM),BKROT)
                          BK(1)=BKROT(1)*BR1(1,1)+BKROT(2)*BR1(1,2)+BKROT(3)*BR1(1,3)
                          BK(2)=BKROT(1)*BR1(2,1)+BKROT(2)*BR1(2,2)+BKROT(3)*BR1(2,3)
                          BK(3)=BKROT(1)*BR1(3,1)+BKROT(2)*BR1(3,2)+BKROT(3)*BR1(3,3)
                          CALL ROTATE (BK,ROTLOC(1,1,JATOM),BKRLOC)
                          CALL YLM (BKRLOC,LOMAX,YL)
                          arg=bkrot(1)*pos(1,lfirst)+  &
                               bkrot(2)*pos(2,lfirst)+  &
                               bkrot(3)*pos(3,lfirst)+  &
                               bkx(i)*transij(1,latom)+bky(i)*transij(2,latom)+bkz(i)*transij(3,latom)
                          phshel=exp((0,1)*(arg*TAU))
                          DO NUM=1,NB
                             PHS(NUM)=PHSHEL*A(I,NUM,kkk)
                          ENDDO
                          DO M=-l,l
                             index=l*(l+1)+m+1
                             DO NUM=1,NB
                                DO irf=1,nrf
                                   ALM(index,num,mu,irf)=ALM(INDEX,num,mu,irf)+ &
                                        ALo(l,jlo,irf)*dconjg(YL(INDEX))*PHS(NUM)
                                ENDDO
                             ENDDO
                          ENDDO
                       ENDDO
                    ENDDO
                 ENDDO
              ENDDO
           ENDIF lo
!*********************************************************
           INDEX=0
           DO L=0,LMAX2
              MAX=2*L+1
              CFAC=(0,1)**L
              DO M=1,MAX
                 INDEX=INDEX+1
                 DO NUM=1,NB
                    DO irf=1,nrf
                       ALM(INDEX,NUM,mu,irf)=ALM(INDEX,NUM,mu,irf)*FAC*CFAC
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
        ENDDO equiv1

        equiv2: DO mu=1,MULT(JATOM)
           latom=lfirst-1+mu
           DO num=1,NB
              INDEX=0
              DO L=0,LMAX2
                 pi12lo=0.d0
                 if (l.le.lomax) then
                    if (loor(l)) CALL RINT13(rf1(1,l,1),rf2(1,l,1), &
                         rf1(1,l,3),rf2(1,l,3),pi12lo,JATOM)
                 endif
                 MAX=2*L+1
                 DO M=1,MAX
                    INDEX=INDEX+1
                    DO NUM2=1,NPROJ
                       projection(num,num2,kkk)=projection(num,num2,kkk) + &
                            conjg(ALM(INDEX,num,mu,1))*C(NUM2,INDEX,LATOM) + &
                            conjg(ALM(INDEX,num,mu,3))*C(NUM2,INDEX,LATOM) * &
                            pi12lo
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
        ENDDO equiv2
     ENDDO kpoints
  ENDDO atoms

  DO KKK=1,NUM_KPTS
     DO num2=1,NPROJ
        DO num1=1,NB
           write(unit_amn, '(2i4,1x,i5,1x,2e18.5)') &
                num1,num2, kkk, projection(num1,num2,KKK) !real(projection(num1,num2,KKK), R8), &
!                aimag(projection(num1,num2,KKK))
        ENDDO
     ENDDO
  ENDDO

  RETURN
end subroutine l2amn
end module     l2amn_m


!!/---
!! Local Variables:
!! mode: f90
!! End:
!!\---
!!
!! Time-stamp: <2016-07-18 15:56:42 assman@faepop71.tu-graz.ac.at>
