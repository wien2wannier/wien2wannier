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
  private; save
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

    allocate( a(Nmat, Nb, num_kpts), vecsz(num_kpts) )

    allocate( BKX(Nmat),         BKY(Nmat),         BKZ(Nmat),          &
         &     XK(     num_kpts), YK(     num_kpts), ZK(     num_kpts), &
         &     GX(Nmat,num_kpts), GY(Nmat,num_kpts), GZ(Nmat,num_kpts))

    GX = 0; GY = 0; GZ = 0
  end subroutine init_xa3
end module xa3


!---------------------------  Procedure modules  ---------------------------
module     almgen_m; contains
subroutine almgen(Alm, stru, jatom, lfirst, NB, kkk)
  use param,     only: lomax, NRF, Lmax2
  use w2w,       only: Nmat, iBlock
  use xa3,       only: vecsz, XK,YK,ZK, BKX,BKY,BKZ, GX,GY,GZ, A
  use xa,        only: fj, dfj, bk, phs
  use lolog,     only: nlo,nlov,nlon,ilo,lapw,n_rad
  use atspdt,    only: P, DP
  use loabc,     only: alo
  use const,     only: R8, C16, TAU
  use structmod, only: struct_t

  !! procedure includes
  use harmon_m
  use Ylm_m

  implicit none

  type(struct_t), intent(in)  :: stru
  integer,        intent(in)  :: jatom, lfirst, NB, kkk
  complex(C16),   intent(out) :: Alm(Nb, Nrf, (LMAX2+1)**2, stru%mult(jatom))

  real(R8) :: h_al(iblock), h_bl(iblock)

  complex(C16) :: PHSHEL,CFAC
  complex(C16) :: YL    ((LMAX2+1)**2)
  complex(C16) :: h_yl  ((LMAX2+1)**2, iblock)
  complex(C16) :: h_alyl((LMAX2+1)**2, iblock)
  complex(C16) :: h_blyl((LMAX2+1)**2, iblock)
  complex(C16) :: ALMt  ((LMAX2+1)**2, Nb, Nrf, stru%mult(JATOM))

  real(R8) :: fac, arg, arg1, arg2, arg3
  integer  :: N, num, i,ii,i3, latom, mu, index, l, max, m,m1, irf,mirf, ibb, lda,ldb,ldc, jlo,jneq

  FAC=2*TAU*stru%RMT(JATOM)**2/sqrt(stru%vol)

  N=vecsz(kkk)
  ALMt=0
  do I=1,N
     BKX(I)=XK(kkk)+GX(I,kkk)                    ! plane wave Q+k
     BKY(I)=YK(kkk)+GY(I,kkk)
     BKZ(I)=ZK(kkk)+GZ(I,kkk)
  end do

  ! calculate Bessel functions j_l(|Q+k|*RMT)
  call harmon(stru, jatom, N, BKX, BKY, BKZ ,Lmax2, FJ, DFJ)
  do MU=1, stru%mult(jatom)
     latom=lfirst-1+mu

     ! matching of plane wave and atomic functions
     do ii=1,N-(nlo+nlon+nlov),iblock
        i3=0
        do i=ii,min(ii+iblock-1,N-(nlo+nlon+nlov))
           i3=i3+1

           ! transform Q+k to cartesian coordinates
           BK = matmul(stru%gbas, (/BKX(I), BKY(I), BKZ(I)/))

           call YLM (BK,LMAX2,YL) ! calculate Y_lm(Q+k)

           ARG1=BKX(I) * stru%POS(1,LATOM)*TAU
           ARG2=BKY(I) * stru%POS(2,LATOM)*TAU
           ARG3=BKZ(I) * stru%POS(3,LATOM)*TAU

           PHSHEL=exp((0,1)*(ARG1+ARG2+ARG3)) ! exp(i(Q+k)R_j)
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
                 h_AL(i3)=FJ(L,I)/P(L,1)/stru%RMT(JATOM)**2
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

        ibb = min(iblock, N - (nlo+nlon+nlov) - ii + 1)
#ifndef _COMPLEX_
        lda=2*(LMAX2+1)**2
#else
        lda=(LMAX2+1)**2
#endif
        ldc=lda
        ldb=nmat
#ifndef _COMPLEX_
        call dgemm('N', 'N', 2*index, Nb, ibb,           &
             &     1._R8, h_alyl, lda, a(ii,1,kkk), ldb, &
             &     1._R8, almt(1,1,1,mu), ldc)
        call dgemm('N', 'N', 2*index, Nb, ibb,           &
             &     1._R8, h_blyl, lda, a(ii,1,kkk), ldb, &
             &     1._R8, almt(1,1,2,mu), ldc)
#else
        call zgemm('N', 'N', index, Nb, ibb,                  &
             &     (1._R8, 0), h_alyl, lda, a(ii,1,kkk), ldb, &
             &     (1._R8, 0), almt(1,1,1,mu), ldc)
        call zgemm('N', 'N', index, Nb, ibb,                  &
             &     (1._R8, 0), h_blyl, lda, a(ii,1,kkk), ldb, &
             &     (1._R8, 0), almt(1,1,2,mu), ldc)
#endif
     end do

     if (nlo.ne.0) then
        i=n-(nlo+nlon)
        do L=0,LoMAX
           mirf=n_rad(l)
           do jlo=1,ilo(l)
              do jneq=1,stru%mult(jatom)
                 do M1=-l,l
                    i=i+1
                    BK = matmul(stru%gbas, (/BKX(I), BKY(I), BKZ(I)/))

                    call YLM (BK,LOMAX,YL)
                    ARG = dot_product((/BKX(I), BKY(I), BKZ(I)/), &
                         &            stru%pos(:,latom))
                    PHSHEL=exp((0,1)*(ARG*TAU))
                    do NUM=1,NB
                       PHS(NUM)=PHSHEL*A(I,NUM,kkk)
                    enddo
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
subroutine l2Mmn(stru, inwf, num_kpts, NNTOT)
  use const,     only: R8, C16, TAU
  use param,     only: unit_out, unit_vsp, unit_vector, NRF, Lmax2
  use xa,        only: bk
  use xa3,       only: xk, yk, zk
  use bessel,    only: ri_mat
  use Amn_Mmn,   only: overlap
  use pairs,     only: kp, kpb, bqx, bqy, bqz
  use lolog,     only: n_rad
  use structmod, only: struct_t
  use inwfmod,   only: inwf_t

  !! procedure includes
  use almgen_m
  use gaunt1_m
  use radint_m
  use atpar_m
  use Ylm_m

  implicit none

  type(struct_t), intent(in) :: stru
  type(inwf_t),   intent(in) :: inwf
  integer,        intent(in) :: num_kpts, NNtot

  integer  :: pair_index, r_index, iscf, jatom, lfirst, itap, itape, jtape
  integer  :: k1, k2, kkk, Nb
  real(R8) :: tALM, tMeas1, tMeas2, t1, tt0, tt1, fac, bx, by, bz
  real(R8) :: KX1,KY1,KZ1     ! k-points in u.c. coordinates k and k+b

  complex(C16) :: YLB((LMAX2+1)**2) ! spherical harmonics expansion of b
  complex(C16) :: PHSHEL, tmp, tmp1

  complex(C16), allocatable :: alm(:,:,:,:),blm(:,:,:,:)

  real(R8) :: BM, tt2,tt3, t2, arg1, arg2, arg3
  integer  :: index1,index2,indexj, irf1,irf2, lj, m1,m2,mj, mu
  integer  :: latom, l_index, l1, l2

  integer :: k1_prog_itvl
  k1_prog_itvl = min(max(num_kpts/10, 1), 100)

  Nb = inwf%bmax - inwf%bmin + 1

  read(unit_vsp, '(50X,I2,//)') ISCF
  LFIRST=1
  atoms: do jatom = 1, stru%Nneq
     write(unit_out, "(/, '===== atom', I5, ' /', I5, ' =====' /)") &
          jatom, stru%Nneq

     allocate(Alm(Nb, Nrf, (Lmax2+1)**2, stru%mult(jatom)),  &
          &   Blm(Nb, Nrf, (Lmax2+1)**2, stru%mult(jatom)))

     talm=0.
     tmeas1=0.
     tmeas2=0.
     call cputim(t1)

     if(jatom > 1) lfirst = lfirst + stru%mult(jatom-1)
     ITAP=30+JATOM
     itape=unit_vector
     jtape=unit_vsp
     rewind(itape)
     call atpar(stru, JATOM, itape, jtape)
     fac=2*TAU * stru%RMT(jatom)**2 / sqrt(stru%vol)
     rewind(itape)

     pair_index=0
     k1loop: do k1=1,num_kpts      !   loop over k1
        pair_index=pair_index+1
        kkk=KP(pair_index)
        KX1=XK(kkk)
        KY1=YK(kkk)
        KZ1=ZK(kkk)

        call almgen(Alm, stru, jatom, lfirst, Nb, kkk)

        k2loop: do k2=1,NNTOT
           if (k2 > 1) pair_index=pair_index+1
           kkk=KPB(pair_index)
           call cputim(tt0)
           call almgen(Blm, stru, jatom, lfirst, Nb, kkk)
           BLM=conjg(BLM)

           call cputim(tt1)

           ! calculate b=k2-k1 add BQ if going around BZ
           BX=XK(kkk)-KX1+BQX(pair_index)
           BY=YK(kkk)-KY1+BQY(pair_index)
           BZ=ZK(kkk)-KZ1+BQZ(pair_index)

           ! transform to cartesian coordinates
           BK = matmul(stru%gbas, (/BX, BY, BZ/))
           BM = sqrt(sum(BK**2))

           ! compute radial intergrals <R(r)|j_(|b|*r)|R'(r)>
           call radint(stru, jatom, inwf%LJmax, bm)
           ! computer Y_lm(b)
           call YLM (BK, inwf%LJmax, YLB)

           indexj=0
           do LJ = 0, inwf%LJmax
              do MJ=-LJ,LJ
                 indexj=indexj+1
                 YLB(indexj)=conjg(YLB(indexj))*(0,1)**LJ
              enddo
           enddo
           call cputim(tt2)
           muloop: do mu=1, stru%mult(jatom)
              latom=lfirst-1+mu
              arg1=Bx * stru%pos(1,latom)*TAU
              arg2=By * stru%pos(2,latom)*TAU
              arg3=Bz * stru%pos(3,latom)*TAU
              PHSHEL=exp((0,1)*(ARG1+ARG2+ARG3))*2*TAU

              L_index=0
              l1loop: do L1=0,LMAX2
                 l2loop: do L2=0,LMAX2
                    ljloop: do LJ = abs(L1-L2), min(L1+L2, inwf%LJmax), 2
                       if (mod((L1+L2+LJ),2) == 1) cycle
                       if ((L1+L2-LJ) < 0         &
                            .or. (L1-L2+LJ) < 0   &
                            .or. (-L1+L2+LJ) < 0) &
                            & cycle
                       L_index=L_index+1
                       m1loop: do M1=-L1,L1
                          mjloop: do MJ = max(-LJ, -L2-M1), min(LJ, L2-M1)
                             M2=M1+MJ                   ! abs(m2) <= l2 !
                             index1=L1*(L1+1)+M1+1
                             index2=L2*(L2+1)+M2+1
                             indexj=LJ*(LJ+1)+MJ+1
                             tmp=YLB(indexj)*PHSHEL*GAUNT1(L2,LJ,L1,M2,MJ,M1)
                             R_index=0
                             do irf1=1,n_rad(L1)
                                do irf2=1,n_rad(L2)
                                   R_index=R_index+1
                                   tmp1=ri_mat(R_index,L_index)*tmp
                                   call zgeru(NB,NB,tmp1,BLM(1,irf2,index2,mu),1,ALM(1,irf1,index1,mu),1,overlap(1,1,pair_index),size(overlap,1))
                                enddo
                             enddo
                          enddo mjloop
                       enddo m1loop
                    enddo ljloop
                 enddo l2loop
              enddo l1loop
           enddo muloop
           call cputim(tt3)
           talm=talm+tt1-tt0
           tmeas1=tmeas1+tt2-tt1
           tmeas2=tmeas2+tt3-tt2
        end do k2loop

        if (mod(k1, k1_prog_itvl) == 0) &
             write(unit_out, "('k1=', I5, ' /', I5, ' (', I3, '%)')") &
             &    k1, num_kpts, (100*k1)/num_kpts
     end do k1loop
     deallocate(ALM,BLM)

     call cputim(t2)
     talm=talm*(NNTOT+1)/NNTOT
     write(unit_out,*) 'CPU time used for atom ',JATOM,' ='&
          ,t2-t1,talm,tmeas1,tmeas2

        if (mod(k1, k1_prog_itvl) == 0) &
             write(unit_out, "('k1=', I5, ' /', I5, ' (', I3, '%)')") &
             &    k1, num_kpts, (100*k1)/num_kpts
  end do atoms
end subroutine l2MMN
end module     l2Mmn_m


module     l2Amn_m; contains
subroutine l2Amn(stru, inwf, num_kpts)
  use param,     only: Lmax2, NRF, lomax, unit_vsp, unit_out, unit_vector
  use w2w,       only: Nmat, iblock, unit_amn
  use const,     only: R8, C16, TAU
  use xa,        only: BK, FJ, DFJ, PHS, BKrot, BKRloc
  use xa3,       only: A, XK,YK,ZK, GX,GY,GZ, BKX,BKY,BKZ, vecsz
  use lolog,     only: nlo,nlov,nlon,loor,ilo,lapw
  use loabc,     only: alo
  use atspdt,    only: P, DP
  use const,     only: R8, C16
  use radfu,     only: RF1, RF2
  use structmod, only: struct_t
  use inwfmod,   only: inwf_t

  !! procedure includes
  use rint13_m
  use harmon_m
  use atpar_m
  use Ylm_m

  implicit none

  type(struct_t), intent(in) :: stru
  type(inwf_t),   intent(in) :: inwf
  integer,        intent(in) :: num_kpts

  complex(C16) :: Yl((Lmax2+1)**2)        ! list of spherical harmics
  complex(C16) :: projection(inwf%bmax-inwf%bmin+1, inwf%Nproj, num_kpts)
  complex(C16) :: phshel,cfac

  real(R8)     :: h_al(iblock), h_bl(iblock)
  complex(C16) :: h_yl  ( (Lmax2+1)**2, iblock)
  complex(C16) :: h_alyl( (Lmax2+1)**2, iblock)
  complex(C16) :: h_blyl( (Lmax2+1)**2, iblock)
  complex(C16) :: Alm   ( (Lmax2+1)**2, inwf%bmax-inwf%bmin+1, stru%Nat, Nrf)

  integer  :: i,ii,i3,l, m,mu,m1, n, ipr, ib, iY
  integer  :: jatom,jlo,jneq, idx, ibb,iscf,irf, lda,ldb,ldc, latom,lfirst
  integer  :: kkk, max, itap,itape,jtape, Nb
  real(R8) :: fac, arg, pi12lo(0:Lmax2)

  Nb = inwf%bmax-inwf%bmin+1
  projection = 0

  rewind(unit_vsp)
  read(unit_vsp, '(50X,I2,//)') ISCF
  LFIRST=1
  atoms: do JATOM = 1, stru%Nneq
     if(jatom > 1) lfirst = lfirst + stru%mult(jatom-1)
     ITAP=30+JATOM
     itape=unit_vector
     jtape=unit_vsp
     rewind(itape)
     call atpar(stru, jatom, itape, jtape)
     rewind(itape)
     write(unit_out, *)
     write(unit_out, '(7X,  "rotij", 10X, "transij")')
     write(unit_out, '(6X, 23("-"))')

     pi12lo = 0
     do l = 0,LOmax
        if (loor(l)) call rint13(stru, jatom, rf1(:,l,1), rf2(:,l,1), &
             &                   rf1(:,l,3), rf2(:,l,3), pi12lo(l))
     end do

     kpoints: do KKK=1,NUM_KPTS
        N=vecsz(kkk)

        ! plane wave Q+k
        BKX(1:N) = XK(kkk) + GX(1:N, kkk)
        BKY(1:N) = YK(kkk) + GY(1:N, kkk)
        BKZ(1:N) = ZK(kkk) + GZ(1:N, kkk)

        ! calculate Bessel functions j_l(|Q+k|*RMT)
        call harmon(stru, jatom, N, Bkx, Bky, Bkz, Lmax2, FJ, DFJ)

        fac = 2*TAU * stru%RMT(jatom)**2 / sqrt(stru%vol)
        Alm = 0

        equiv1: do mu = 1, stru%mult(jatom) !  loop over equivalent atoms
           latom=lfirst-1+mu
           if (kkk == 1) then
              do i=1,3
                 write(unit_out, '(4x,3I3,4x,f12.4)') &
                      stru%rotij(i,:,latom), stru%transij(i,latom)
              enddo
              write(unit_out,*)
           endif

           ! matching of plane wave and atomic functions
           do ii=1,N-(nlo+nlon+nlov),iblock
              i3=0
              do i=ii,min(ii+iblock-1,N-(nlo+nlon+nlov))
                 i3=i3+1

                 ! go to coord. of latom
                 bkrot = matmul(stru%rotij(:,:,latom), &
                      &         (/bkx(i), bky(i), bkz(i)/))

                 ! transform Q+k to cartesian coordinates
                 bk = matmul(stru%gbas, bkrot)

                 ! go to local (rotated) coords
                 bkrloc = matmul(stru%locrot(:,:,jatom), bk)

                 ! calculate Y_lm(Q+k)
                 call YLM (bkrloc,LMAX2,YL)
                 arg =  dot_product(bkrot, stru%pos(:, lfirst)) &
                      + dot_product((/bkx(i), bky(i), bkz(i)/), &
                      &             stru%transij(:,latom))
                 phshel=exp((0,1)*arg*TAU)

                 idx=0
                 do l=0,Lmax2
                    do m=1, 2*l+1
                       idx=idx+1
                       h_yl(idx,i3)=conjg(yl(idx))*phshel
                    enddo
                 enddo
              enddo

              IDX=0
              do L=0,LMAX2
                 i3=0
                 do i=ii,min(ii+iblock-1,N-(nlo+nlon+nlov))
                    i3=i3+1
                    if (lapw(l)) then
                       h_AL(i3)=DFJ(L,I)*P(L,2)-FJ(L,I)*DP(L,2)
                       h_BL(i3)=FJ(L,I)*DP(L,1)-DFJ(L,I)*P(L,1)
                    else
                       h_AL(i3)=FJ(L,I)/P(L,1)/stru%RMT(jatom)**2
                       h_BL(i3) = 0
                    endif
                 enddo
                 MAX=2*L+1
                 do M=1,MAX
                    IDX=IDX+1
                    i3=0
                    do i=ii,min(ii+iblock-1,N-(nlo+nlon+nlov))
                       i3=i3+1
                       h_alyl(idx,i3)=h_AL(i3)*h_YL(idx,i3)
                       h_blyl(idx,i3)=h_BL(i3)*h_YL(idx,i3)
                    enddo
                 enddo
              enddo
              ibb=min(iblock,N-(nlo+nlon+nlov)-ii+1)
#ifndef _COMPLEX_
              lda=2*(LMAX2+1)*(LMAX2+1)
#else
              lda=(LMAX2+1)*(LMAX2+1)
#endif
              ldc=lda
              ldb=nmat

#ifndef _COMPLEX_
              call dgemm('N','N',2*idx,nb,ibb,1.d0, &
                   h_alyl,lda,a(ii,1,kkk),ldb,1.d0, &
                   alm(1,1,mu,1),ldc)
              call dgemm('N','N',2*idx,nb,ibb,1.d0, &
                   h_blyl,lda,a(ii,1,kkk),ldb,1.d0, &
                   alm(1,1,mu,2),ldc)
#else
              call zgemm('N','N',idx,nb,ibb,(1.d0,0.d0), &
                   h_alyl,lda,a(ii,1,kkk),ldb,(1.d0,0.d0), &
                   alm(1,1,mu,1),ldc)
              call zgemm('N','N',idx,nb,ibb,(1.d0,0.d0), &
                   h_blyl,lda,a(ii,1,kkk),ldb,(1.d0,0.d0), &
                   alm(1,1,mu,2),ldc)
#endif
           enddo
!********  LOCAL ORBITALS *****************************
           lo: if (nlo.ne.0) then
              i=n-(nlo+nlon)
              do L=0,LoMAX
                 do jlo=1,ilo(l)
                    do jneq=1,stru%mult(jatom)
                       do M1=-l,l
                          i=i+1
                          BK(1)=BKX(I)
                          BK(2)=BKY(I)
                          BK(3)=BKZ(I)
                          bkrot = matmul(stru%rotij(:,:,latom), BK)
                          BK    = matmul(stru%gbas, BKROT)

                          bkrloc = matmul(stru%locrot(:,:,jatom), BK)
                          call Ylm (Bkrloc, LOmax, YL)
                          arg =  dot_product(bkrot(:), stru%pos(:, lfirst)) &
                               + dot_product((/bkx(i), bky(i), bkz(i)/),    &
                               &             stru%transij(:,latom))

                          phshel=exp((0,1)*TAU*arg)
                          PHS(:) = phshel * A(i,:,kkk)

                          do m=-l,l
                             idx=l*(l+1)+m+1
                             do ib=1,Nb
                                Alm(idx, ib, mu, :)=Alm(idx, ib, mu, :) + &
                                     Alo(l, jlo, :) * conjg(Yl(idx))*phs(ib)
                             end do
                          end do
                       end do
                    end do
                 end do
              end do
           end if lo
!*********************************************************
           idx=0
           do l=0,Lmax2
              cfac = (0,1)**l
              do m=1, 2*l+1
                 idx=idx+1
                 do ib=1,Nb
                    do irf=1,Nrf
                       Alm(idx,ib,mu,irf)=Alm(idx,ib,mu,irf)*fac*cfac
                    enddo
                 enddo
              enddo
           enddo
        enddo equiv1

        equiv2: do mu=1,stru%mult(jatom)
           latom = lfirst - 1 + mu

           do ib = 1,Nb
              do ipr = 1, inwf%Nproj
                 do iY = 1, inwf%projections(ipr)%NY
                    if (inwf%projections(ipr)%iat(iY) /= latom) &
                         cycle
                    l = inwf%projections(ipr)%l(iY)
                    m = inwf%projections(ipr)%m(iY)
                    ! idx is the compound index to Alm
                    idx = l**2 + m + l + 1

                    projection(ib, ipr, kkk) = projection(ib, ipr, kkk) &
                         & + conjg(Alm(idx, ib, mu, 1))                 &
                         &   * inwf%projections(ipr)%coeff(iY)          &
                         & + conjg(Alm(idx, ib, mu, 3))                 &
                         &   * inwf%projections(ipr)%coeff(iY)          &
                         &   * pi12lo(l)
                 enddo
              enddo
           enddo
        enddo equiv2
     enddo kpoints
  enddo atoms

  do kkk = 1,num_kpts
     do ipr = 1, inwf%Nproj
        do ib = 1,Nb
           write(unit_amn, '(2i4,1x,i5,1x,2e18.5)') &
                ib, ipr, kkk, projection(ib, ipr, KKK)
        enddo
     enddo
  enddo

  return
end subroutine l2Amn
end module     l2Amn_m


module     planew_m; contains
subroutine planew(stru, NB, num_kpts, Nntot, maxx, maxy, maxz)
  use param,     only: unit_out
  use w2w,       only: unit_mmn
  use xa3,       only: gx,gy,gz, a, vecsz
  use amn_mmn,   only: overlap
  use pairs,     only: kp, kpb, bqx,bqy,bqz, bqx1,bqy1,bqz1
  use lolog,     only: nlo,nlov,nlon
  use const,     only: R8, C16, TAU
  use structmod, only: struct_t

  implicit none

  type(struct_t), intent(in) :: stru
  integer,        intent(in) :: NB, num_kpts, NNTOT, maxx,maxy,maxz

  integer  :: index, pair_index, m2x,m2y,m2z, max, i,j,k, ii,jj,kk
  integer  :: i_x,i_y,i_z, k1,k2, kkk, lfirst, jatom, latom, mu
  integer  :: num1,num2, id, ix1,ix2,iy1,iy2,iz1,iz2
  real(R8) :: t1,t2, tt1,tt2,tt3, t_b, t_prod, vv, bx,by,bz, x, arg

  logical, allocatable, dimension(:,:) :: a_null, a_null_fin
  complex(C16)       sum(NB),SUMA(NB,NB),FF
  complex(c16), allocatable, dimension(:,:,:) :: FAC
  complex(C16)       PHSHEL
  complex(C16), allocatable, dimension(:,:,:) :: A_
  complex(C16), allocatable, dimension(:,:) :: B_

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

        A_(:,index,kkk)=A(I,:,kkk)
     end do
     do k2=1,NNTOT
        if (k2 > 1) pair_index=pair_index+1
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
  do JATOM=1,stru%Nneq
     if (jatom > 1) lfirst = lfirst + stru%mult(jatom-1)
     VV=2*TAU * stru%RMT(jatom)**3/3
     do K=-m2z,m2z
        do J=-m2y,m2y
           do I=-m2x,m2x
              BX = dot_product(stru%gbas(1, :), (/ i, j, k /))
              BY = dot_product(stru%gbas(2, :), (/ i, j, k /))
              BZ = dot_product(stru%gbas(3, :), (/ i, j, k /))

              X = stru%RMT(JATOM) * sqrt(BX**2 + BY**2 + BZ**2)
              phshel=0
              do mu=1,stru%mult(jatom)
                 latom=lfirst+mu-1
                 arg =  dot_product((/i,j,k/), stru%pos(:,LATOM))
                 phshel = phshel + exp(-(0,1)*arg*TAU)
              end do

              if (i==0 .and. j==0 .and. k==0) then
                 fac(i,j,k) = fac(i,j,k) - stru%mult(jatom) * VV/stru%vol
              else
                 fac(i,j,k) = fac(i,j,k) &
                      - 3*phshel * (sin(x) - x*cos(x)) / x**3 * VV/stru%vol
              end if
           end do
        end do
     end do
  end do
  FAC(0,0,0)=FAC(0,0,0)+1
  call cputim(t2)
  write(unit_out,*)'t_init=',t2-t1

  ! read in k-point and basis vectors
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
        if (k2 > 1) pair_index=pair_index+1
        kkk=KPB(pair_index)

        id=-(BQZ(pair_index)*m2y+BQY(pair_index))*m2x-BQX(pair_index)-1
        SUMA=0

        do J=0,max
           id=id+1

           if (a_null(j,kkk) .and. 0 <= id .and. id < max) then
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
     if (k1 == 1) write(unit_out,*)'PLANEW TIME PER k-point:',t_b+t_prod
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
subroutine read_vec(NEmin, NEmax, Nneq, kkk, maxx, maxy, maxz, Efermi)
  use param,  only: unit_out, unit_vector
  use w2w,    only: unit_eig
  use const,  only: R8, Ryd_eV
  use xa3,    only: XK,YK,ZK, GX,GY,GZ, vecsz, A

  implicit none

  integer, intent(in)    :: NEmin, NEmax, Nneq
  integer, intent(out)   :: maxx, maxy, maxz
  integer, intent(inout) :: kkk
  real(R8),intent(in)    :: efermi

  integer :: i, j, n, NB, NE, num

  real(r8) :: E(1000)

  character(10)    BNAME

  do I = 1, Nneq
     read(unit_vector)
     read(unit_vector)
  enddo

  maxx=0; maxy=0; maxz=0
  kpoint: do
     read(unit_vector,end=998) XK(kkk+1),YK(kkk+1),ZK(kkk+1),BNAME,N,NE
     KKK=KKK+1
     vecsz(kkk)=N
     read(unit_vector) (GX(I,kkk),GY(I,kkk),GZ(I,kkk), I=1,N)
     do I=1,N
        if (abs(GX(I,kkk)) > maxx) maxx=abs(GX(I,kkk))
        if (abs(GY(I,kkk)) > maxy) maxy=abs(GY(I,kkk))
        if (abs(GZ(I,kkk)) > maxz) maxz=abs(GZ(I,kkk))
     enddo

     do J=1,NE
        read(unit_vector) num, E(num)
        if (num >= NEmin .and. num <= NEmax) then
           read(unit_vector)(A(I,NUM-NEMIN+1,kkk),I=1,N)
        else
           read(unit_vector)
        endif
     enddo

     do NUM=NEMIN,NEMAX
        NB=NUM-NEMIN+1
        write(unit_eig,"(2I12,F22.16)")NB,kkk, Ryd_eV*(E(NUM)-efermi)
     enddo
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
!! Time-stamp: <2016-08-04 09:19:21 assman@faepop71.tu-graz.ac.at>
