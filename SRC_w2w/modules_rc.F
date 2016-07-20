!!! wien2wannier/SRC_w2w/modules_rc.F
!!!
!!!    Modules for wien2wannier.  This file contains modules that
!!!    depend on real/complex compilation.
!!!
!!! Copyright 2010-2012 Jan Kuneš, Philipp Wissgott
!!!           2013-2016 Elias Assmann

!!/=== Defined here: =============================
!!
!! xa3: init_xa3(), BKX(:),BKY(:),BKZ(:), GX(:,:),GY(:,:),GZ(:,:),
!!      XK(:),YK(:),ZK(:), vecsz(:)
!!
!! Procedure modules:
!!
!!    almgen_m, l2Mmn_m, l2Amn_m, planew_m, read_vec_m
!!
!!\===============================================

module xa3
  use const, only: R8
#ifdef _COMPLEX_
  use const, only: C16
#endif

  implicit none
  private
  public :: init_xa3, XK,YK,ZK, BKX,BKY,BKZ, GX,GY,GZ, vecsz, a

#ifndef _COMPLEX_
  real(R8), &
#else
  complex(C16), &
#endif
  allocatable :: a(:,:,:)

  integer,  allocatable :: GX(:,:),GY(:,:),GZ(:,:), vecsz(:)
  real(R8), allocatable :: BKX(:), BKY(:), BKZ(:)
  real(R8), allocatable :: XK(:),  YK(:),  ZK(:)

contains
  subroutine init_xa3(nb, nmat, num_kpts)
    integer, intent(in) :: NB, Nmat, num_kpts

    allocate( a  (Nmat, Nb, num_kpts), vecsz(num_kpts) )
    allocate( BKX(Nmat),         BKY(Nmat),         BKZ(Nmat),          &
         &    GX (Nmat,num_kpts), GY(Nmat,num_kpts), GZ(Nmat,num_kpts), &
         &    XK (     num_kpts), YK(     num_kpts), ZK(     num_kpts))

    GX = 0; GY = 0; GZ = 0
  end subroutine init_xa3
end module xa3


!---------------------------  Procedure modules  ---------------------------
module     almgen_m; contains
subroutine almgen(Alm,jatom,lfirst,NB,kkk)
  use param,  only: lomax, NRF, Lmax2
  use w2w,    only: Nmat, iBlock
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


module     l2Mmn_m; contains
subroutine l2Mmn(NB,num_kpts,NNTOT,LJMAX)
  use const,   only: R8, C16, TAU
  use param,   only: unit_out, unit_vsp, unit_vector, NRF, Lmax2
  use struct,  only: Nat, mult, RMT, vol, pos
  use xa,      only: bk
  use xa3,     only: xk, yk, zk
  use bessel,  only: ri_mat
  use Amn_Mmn, only: overlap
  use pairs,   only: kp, kpb, bqx, bqy, bqz
  use gener,   only: br1
  use lolog,   only: n_rad

  !! procedure includes
  use almgen_m
  use gaunt1_m
  use radint_m
  use atpar_m
  use Ylm_m

  implicit none

  integer, intent(in) :: Nb, num_kpts, NNtot, LJmax

  integer  :: pair_index, r_index, iscf, jatom, lfirst, itap, itape, jtape
  integer  :: k1, k2, kkk
  real(R8) :: tALM, tMeas1, tMeas2, t1, tt0, tt1, fac, bx, by, bz
  real(R8) :: KX1,KY1,KZ1     ! k-points in u.c. coordinates k and k+b

  COMPLEX(C16)       YLB((LMAX2+1)*(LMAX2+1))       ! spherical harmonics expansion of b
  COMPLEX(C16)       PHSHEL,tmp,tmp1

  !...............................
  complex(C16), allocatable :: alm(:,:,:,:),blm(:,:,:,:)
  !...............................

  real(R8) :: BM, tt2,tt3, t2, arg1, arg2, arg3
  integer  :: index1,index2,indexj, irf1,irf2, lj, m1,m2,mj, mu, latom, l_index, l1, l2

  integer :: k1_prog_itvl
  k1_prog_itvl = min(max(num_kpts/10, 1), 100)

  !------------------------------------------------------------------

  !     ------------------
  !     LOOP FOR ALL ATOMS
  !     ------------------

  READ(unit_vsp,2032) ISCF
  LFIRST=1
  atoms: DO JATOM=1,NAT
     write(unit_out, "(/, '===== atom', I5, ' /', I5, ' =====' /)") jatom, nat

     ALLOCATE(ALM(NB,NRF,(LMAX2+1)*(LMAX2+1),MULT(JATOM)),  &
          BLM(NB,NRF,(LMAX2+1)*(LMAX2+1),MULT(JATOM)))
     talm=0.
     tmeas1=0.
     tmeas2=0.
     call cputim(t1)
     IF(JATOM.GT.1) LFIRST=LFIRST+MULT(JATOM-1)
     ITAP=30+JATOM
     itape=unit_vector
     jtape=unit_vsp
     rewind(itape)
     call atpar(JATOM, itape, jtape)
     FAC=2*TAU*RMT(JATOM)**2/sqrt(VOL)
     rewind(itape)

     !................................
     pair_index=0
     k1loop: DO k1=1,num_kpts      !   loop over k1
        pair_index=pair_index+1
        kkk=KP(pair_index)
        KX1=XK(kkk)
        KY1=YK(kkk)
        KZ1=ZK(kkk)
        !	write(92,*)'k-point:',kkk,KX1,KY1,KZ1
        call almgen(ALM,JATOM,LFIRST,NB,KKK)
        k2loop: DO k2=1,NNTOT
           if (k2.gt.1) pair_index=pair_index+1
           kkk=KPB(pair_index)
           call cputim(tt0)
           CALL almgen(BLM,JATOM,LFIRST,NB,KKK)
           BLM=conjg(BLM)
           call cputim(tt1)
           BX=XK(kkk)-KX1+BQX(pair_index)                 ! calculate b=k2-k1 add BQ if going around BZ
           BY=YK(kkk)-KY1+BQY(pair_index)
           BZ=ZK(kkk)-KZ1+BQZ(pair_index)
           !	 write(92,*)BX,BY,BZ
           BK(1)=BX*BR1(1,1)+BY*BR1(1,2)+BZ*BR1(1,3)      ! transform to cartesian coordinates
           BK(2)=BX*BR1(2,1)+BY*BR1(2,2)+BZ*BR1(2,3)
           BK(3)=BX*BR1(3,1)+BY*BR1(3,2)+BZ*BR1(3,3)
           BM=SQRT(BK(1)*BK(1)+BK(2)*BK(2)+BK(3)*BK(3))
           call radint(JATOM,LJMAX,BM)                  ! compute radial intergrals <R(r)|j_(|b|*r)|R'(r)>
           CALL YLM (BK,LJMAX,YLB)                      ! computer Y_lm(b)
           indexj=0
           DO LJ=0,LJMAX
              DO MJ=-LJ,LJ
                 indexj=indexj+1
                 YLB(indexj)=conjg(YLB(indexj))*(0,1)**LJ
              ENDDO
           ENDDO
           call cputim(tt2)
           muloop: DO mu=1,MULT(JATOM)
              latom=lfirst-1+mu
              ARG1=BX*POS(1,LATOM)*TAU
              ARG2=BY*POS(2,LATOM)*TAU
              ARG3=BZ*POS(3,LATOM)*TAU
              PHSHEL=exp((0,1)*(ARG1+ARG2+ARG3))*2*TAU

              ! overlap=conjg(alm(2))*alm(1)*gaunt(2,j,1)*rad_int(2,j,1)
              L_index=0
              l1loop: DO L1=0,LMAX2
                 l2loop: do L2=0,LMAX2
                    ljloop: DO LJ = abs(L1-L2), min(L1+L2, LJMAX), 2
           IF (MOD((L1+L2+LJ),2) .EQ. 1) cycle
           IF ((L1+L2-LJ).LT.0.OR.(L1-L2+LJ).LT.0.OR.(-L1+L2+LJ).LT.0) cycle
                       L_index=L_index+1
                       m1loop: DO M1=-L1,L1
                          mjloop: DO MJ = max(-LJ, -L2-M1), min(LJ, L2-M1)
                             M2=M1+MJ                   ! abs(m2) <= l2 !
                             index1=L1*(L1+1)+M1+1
                             index2=L2*(L2+1)+M2+1
                             indexj=LJ*(LJ+1)+MJ+1
                             tmp=YLB(indexj)*PHSHEL*GAUNT1(L2,LJ,L1,M2,MJ,M1)
                             R_index=0
                             DO irf1=1,n_rad(L1)
                                do irf2=1,n_rad(L2)
                                   R_index=R_index+1
                                   tmp1=ri_mat(R_index,L_index)*tmp
                                   call zgeru(NB,NB,tmp1,BLM(1,irf2,index2,mu),1,ALM(1,irf1,index1,mu),1,overlap(1,1,pair_index),size(overlap,1))
                                ENDDO
                             ENDDO
                          ENDDO mjloop
                       ENDDO m1loop
                    ENDDO ljloop
                 ENDDO l2loop
              ENDDO l1loop
           ENDDO muloop
           call cputim(tt3)
           talm=talm+tt1-tt0
           tmeas1=tmeas1+tt2-tt1
           tmeas2=tmeas2+tt3-tt2
        END DO k2loop

        if (mod(k1, k1_prog_itvl) == 0) &
             write(unit_out, "('k1=', I5, ' /', I5, ' (', I3, '%)')") &
             &    k1, num_kpts, (100*k1)/num_kpts
     END DO k1loop
     DEALLOCATE(ALM,BLM)
     call cputim(t2)
     talm=talm*(NNTOT+1)/NNTOT
     write(unit_out,*) 'CPU time used for atom ',JATOM,' ='&
          ,t2-t1,talm,tmeas1,tmeas2

        if (mod(k1, k1_prog_itvl) == 0) &
             write(unit_out, "('k1=', I5, ' /', I5, ' (', I3, '%)')") &
             &    k1, num_kpts, (100*k1)/num_kpts
  END DO atoms

  ! ....END LOOP OVER ALL ATOMS

2032 format(50X,I2,//)
end subroutine l2MMN
end module     l2Mmn_m


module     l2Amn_m; contains
subroutine l2Amn(NB,NPROJ,NUM_KPTS)
  use param,   only: Lmax2, NRF, lomax, unit_vsp, unit_out, unit_vector
  use w2w,     only: Ndif, Nmat, iblock, unit_amn
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
     itape=unit_vector
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
end subroutine l2Amn
end module     l2Amn_m


module     planew_m; contains
subroutine planew(NB,num_kpts,NNTOT,maxx,maxy,maxz)
  use param,   only: unit_out
  use w2w,     only: unit_mmn
  use struct,  only: Nat, mult, RMT, pos, vol
  use xa3,     only: gx,gy,gz, a, vecsz
  use amn_mmn, only: overlap
  use pairs,   only: kp, kpb, bqx,bqy,bqz, bqx1,bqy1,bqz1
  use lolog,   only: nlo,nlov,nlon
  use gener,   only: br1
  use const,   only: R8, C16, TAU

  implicit none

  integer, intent(in) :: NB, num_kpts, NNTOT, maxx,maxy,maxz

  integer  :: index, pair_index, m2x,m2y,m2z, max, i,j,k, ii,jj,kk
  integer  :: i_x,i_y,i_z, k1,k2, kkk, lfirst, jatom, latom, mu
  integer  :: num,num1,num2, id, ix1,ix2,iy1,iy2,iz1,iz2
  real(R8) :: t1,t2, tt1,tt2,tt3, t_b, t_prod, vv, bx,by,bz, x, arg

  logical, allocatable, dimension(:,:) :: a_null, a_null_fin
  COMPLEX(C16)       SUM(NB),SUMA(NB,NB),FF
  complex(c16), allocatable, dimension(:,:,:) :: FAC
  COMPLEX(C16)       PHSHEL
  COMPLEX(C16), allocatable, dimension(:,:,:) :: A_
  COMPLEX(C16), allocatable, dimension(:,:) :: B_

  !allocate memory to arrays
  allocate(FAC(-(2*maxx+1):(2*maxx+1),-(2*maxy+1):(2*maxy+1),-(2*maxz+1):(2*maxz+1)) )
  allocate(A_(NB,0:(2*maxx+1)*(2*maxy+1)*(2*maxz+1)-1,NUM_KPTS) )
  allocate(B_(NB,0:(2*maxx+1)*(2*maxy+1)*(2*maxz+1)-1) )
  allocate(a_null(0:(2*maxx+1)*(2*maxy+1)*(2*maxz+1)-1,NUM_KPTS) )
  allocate(a_null_fin(0:(2*maxx+1)*(2*maxy+1)*(2*maxz+1)-1,NUM_KPTS) )

  m2x=2*maxx+1
  m2y=2*maxy+1
  m2z=2*maxz+1
  max=m2x*m2y*m2z-1
  write(unit_out,*)'MAX Q:',maxx,maxy,maxz
  t_b=0.d0
  t_prod=0.d0

  call cputim(t1)
  A_=0
  a_null=.false.
  a_null_fin=.false.
  pair_index=0
  do k1=1,num_kpts
     pair_index=pair_index+1
     kkk=KP(pair_index)
     do I=1,vecsz(kkk)-(nlo+nlon+nlov)
        I_X=GX(I,kkk)+maxx
        I_Y=GY(I,kkk)+maxy
        I_Z=GZ(I,kkk)+maxz
        index=(I_Z*m2y+I_Y)*m2x+I_X
        a_null(index,kkk)=.true.
        do num=1,nb
           A_(num,index,kkk)=A(I,num,kkk)
        end do
     end do
     do k2=1,NNTOT
        if (k2.gt.1) pair_index=pair_index+1
        kk=KPB(pair_index)
        do I=1,vecsz(kkk)
           I_X=GX(I,kk)+maxx-BQX(pair_index)
           I_Y=GY(I,kk)+maxy-BQY(pair_index)
           I_Z=GZ(I,kk)+maxz-BQZ(pair_index)
           index=(I_Z*m2y+I_Y)*m2x+I_X
           a_null_fin(index,kkk)=.true.
        end do
     end do
  end do

  FAC=0
  LFIRST=1
  do JATOM=1,NAT
     IF (JATOM.GT.1) LFIRST=LFIRST+MULT(JATOM-1)
     VV=2*TAU*RMT(JATOM)**3/3
     do K=-m2z,m2z
        do J=-m2y,m2y
           do I=-m2x,m2x
              BX=I*BR1(1,1)+J*BR1(1,2)+K*BR1(1,3)
              BY=I*BR1(2,1)+J*BR1(2,2)+K*BR1(2,3)
              BZ=I*BR1(3,1)+J*BR1(3,2)+K*BR1(3,3)
              X=RMT(JATOM)*SQRT(BX*BX+BY*BY+BZ*BZ)
              PHSHEL=0
              do MU=1,MULT(JATOM)
                 LATOM=LFIRST+MU-1
                 ARG=I*POS(1,LATOM)  &
                      +J*POS(2,LATOM)  &
                      +K*POS(3,LATOM)
                 PHSHEL=PHSHEL+exp(-(0,1)*ARG*TAU)
              end do
              if ((I*I+J*J+K*K).eq.0) then
                 FAC(I,J,K)=FAC(I,J,K)-MULT(JATOM)*VV/VOL
              ELSE
                 FAC(I,J,K)=FAC(I,J,K)-3*PHSHEL*(SIN(X)-X*COS(X))/(X*X*X)*VV/VOL
              end if
           end do
        end do
     end do
  end do
  FAC(0,0,0)=FAC(0,0,0)+1.d0
  call cputim(t2)
  write(unit_out,*)'t_init=',t2-t1

  !............................................................................
  !.....READ IN K POINT AND BASIS VECTORS
  pair_index=0
  do k1=1,num_kpts
     pair_index=pair_index+1
     kkk=KP(pair_index)
     call cputim(tt1)
     b_=0
     do J=0,max
        if (a_null_fin(J,kkk)) then
           SUM=0
           JJ=J/m2x
           iz2=jj/m2y
           iy2=jj-iz2*m2y
           ix2=j-jj*m2x
           do I=0,max
              if(a_null(I,kkk)) then
                 II=I/m2x
                 iz1=II/m2y
                 iy1=II-iz1*m2y
                 ix1=I-II*m2x
                 FF=FAC(ix2-ix1,iy2-iy1,iz2-iz1)
                 sum=sum+a_(1:nb,I,k1)*FF
              end if
           end do
           b_(1:nb,j)=sum
        end if
     end do
     call cputim(tt2)
     do k2=1,NNTOT
        if (k2.gt.1) pair_index=pair_index+1
        kkk=KPB(pair_index)
        !bugfix 07/23 p.wissgott
        !    id=-(BQZ(pair_index)*m2y+BQY(pair_index))*m2z-BQX(pair_index)-1
        id=-(BQZ(pair_index)*m2y+BQY(pair_index))*m2x-BQX(pair_index)-1
        SUMA=0
        !debug
        !id=0
        do J=0,max
           id=id+1
           !        id=j-(BQZ(pair_index)*m2y+BQY(pair_index))*m2z-BQX(pair_index)
           if (a_null(j,kkk).and.0.le.id.and.id.lt.max) then
              do NUM1=1,NB
#ifndef _COMPLEX_
                 SUMA(1:nb,num1)=SUMA(1:nb,num1)+a_(1:nb,j,kkk)*b_(num1,id)
#else
                 SUMA(1:nb,num1)=SUMA(1:nb,num1)+conjg(a_(1:nb,j,kkk))*b_(num1,id)
#endif
              end do
           endif
        end do


        overlap(1:nb,1:nb,pair_index) = overlap(1:nb,1:nb,pair_index) + SUMA
        write(unit_mmn,235)KP(pair_index),KPB(pair_index),BQX1(pair_index),BQY1(pair_index), &
             BQZ1(pair_index)

        do NUM2=1,NB
           do NUM1=1,NB
              write(unit_mmn,236) conjg(overlap(num2,num1,pair_index))
           end do
        end do
     end do
     call cputim(tt3)
     t_b=t_b+tt2-tt1
     t_prod=t_prod+tt3-tt2
     if (k1.eq.1) write(unit_out,*)'PLANEW TIME PER k-point:',t_b+t_prod
  end do
  write(unit_out,*)'t_b=',t_b
  write(unit_out,*)'t_prod=',t_prod


235 format(5i8)
236 format(32f18.12)
  return

  !release the memory
  deallocate(FAC)
  deallocate(A_)
  deallocate(B_)
  deallocate(a_null)
  deallocate(a_null_fin)
end subroutine PLANEW
end module planew_m


module     read_vec_m; contains
!!!    Read a ‘vector’ file
!!!
!!! Copyright 2010-2012 Jan Kuneš, Philipp Wissgott
!!!           2013-2015 Elias Assmann
subroutine read_vec(NEMIN,NEMAX,kkk,maxx,maxy,maxz,efermi)
  use param,  only: unit_out, unit_vector
  use w2w,    only: unit_eig
  use const,  only: R8, Ryd_eV
  USE struct, only: nat
  USE xa3,    only: XK,YK,ZK, GX,GY,GZ, vecsz, A

  implicit none

  integer, intent(in)    :: NEMIN, NEMAX
  integer, intent(inout) :: maxx, maxy, maxz
  integer, intent(inout) :: kkk
  real(R8),intent(in)    :: efermi

  integer :: i, j, n, NB, NE, num

  real(r8) :: E(1000)

  CHARACTER(10)    BNAME

  DO I=1,NAT
     READ(unit_vector)
     READ(unit_vector)
  ENDDO

  maxx=0; maxy=0; maxz=0
  kpoint: do
     READ(unit_vector,END=998) XK(kkk+1),YK(kkk+1),ZK(kkk+1),BNAME,N,NE
     KKK=KKK+1
     vecsz(kkk)=N
     READ(unit_vector) (GX(I,kkk),GY(I,kkk),GZ(I,kkk), I=1,N)
     DO I=1,N
        IF (abs(GX(I,kkk)).gt.maxx) maxx=abs(GX(I,kkk))
        IF (abs(GY(I,kkk)).gt.maxy) maxy=abs(GY(I,kkk))
        IF (abs(GZ(I,kkk)).gt.maxz) maxz=abs(GZ(I,kkk))
     ENDDO

     DO J=1,NE
        READ(unit_vector)NUM,E(NUM)
        if (NUM.ge.NEMIN.and.NUM.le.nemax) then
           READ(unit_vector)(A(I,NUM-NEMIN+1,kkk),I=1,N)
        else
           READ(unit_vector)
        endif
     ENDDO

     DO NUM=NEMIN,NEMAX
        NB=NUM-NEMIN+1
        write(unit_eig,"(2I12,F22.16)")NB,kkk, ryd_ev*(E(NUM)-efermi)
     ENDDO
  end do kpoint

998 write(unit_out,*)'vector read in',kkk
end subroutine read_vec
end module read_vec_m

!!/---
!! Local Variables:
!! mode: f90
!! End:
!!\---
!!
!! Time-stamp: <2016-07-20 10:35:44 assman@faepop71.tu-graz.ac.at>
