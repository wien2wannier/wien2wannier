!!! wien2wannier/SRC_w2w/modules.f
!!!
!!!    Modules for wien2wannier.  This file contains modules that are
!!!    independent of real/complex compilation.
!!!
!!! Copyright 2010-2012 Jan Kuneš, Philipp Wissgott
!!!           2013-2016 Elias Assmann

!!/=== Defined here: =============================
!!
!! Amn_Mmn: init_Amn_Mmn(), overlap(:,:,:), c(:,:,:)
!!
!! assleg:  init_assleg(), YR(:,:), N, MAXDIM
!!
!! atspdt:  P(0:Lmax2, Nrf), DP(0:Lmax2,Nrf)
!!
!! bessel:  rj(:,:), ri_mat(:,:), init_bessel()
!!
!! loabc:   alo(0:LOmax, Nloat, Nrf)
!!
!! lolog:   n_rad(0:lmax2), ilo(0:lomax), loor(0:lomax), lapw(0:lmax2),
!!          Nlo, Nlov, Nlon
!!
!! pairs:   init_pairs(), BQX(:), BQY(:), BQZ(:), BQX1(:), BQY1(:), BQZ1(:)
!!          KP(:), KPB(:)
!!
!! radfu:   RF1(Nrad,0:Lmax2,nrf), RF2(Nrad,0:Lmax2,Nrf)
!!
!! w2w:     unit_amn=7,  unit_mmn  =8, unit_nnkp=11, unit_eig=12,
!!          unit_ene=50, unit_fermi=51, iBlock=128, NMAT, NDIF
!!
!! xa:      init_xa(), bk(3), bkrot(3), bkrloc(3), fj(:,:), dfj(:,:),
!!          phs(:), r(:)
!!
!! Procedure modules:
!!
!!    abc_m, atpar_m, gaunt1_m, gaunt2_m, harmon_m, outwin_m,
!!    radint_m, rint13_m
!!
!!\===============================================


module w2w
  implicit none
  public; save

  integer, parameter :: unit_amn=7,  unit_mmn  =8, unit_nnkp=11, unit_eig=12
  integer, parameter :: unit_ene=50, unit_fermi=51
  ! Optimize IBLOCK for your hardware (32-255)
  integer, parameter :: iBlock=128

  integer :: NMAT=0
end module w2w


module assleg
  use const, only: R8

  implicit none
  private; save

  public :: YR, N, init_assleg

  real(R8), allocatable :: YR(:,:)
  integer,  parameter   :: N=6, MAXDIM=81
contains
  subroutine init_assleg
    implicit none

    allocate( YR(N, MAXDIM) )
    YR=0
  end subroutine init_assleg
end module assleg


module bessel
  use const, only: R8
  implicit none
  private :: R8
  public; save

  ! Cave: rj's first index will start at 0
  real(R8), allocatable :: rj(:,:), ri_mat(:,:)

contains
  subroutine init_bessel(LMAX2,LJMAX,NRAD,NRF)
    integer, intent(in) :: Lmax2, LJmax, Nrad, NRF

    integer :: l, l1, l2, lj

    L=0
    do L1=0,LMAX2
       do L2=0,LMAX2
          do LJ=0,LJMAX
             if (mod((L1+L2+LJ),2) .eq. 1) cycle
             if ((L1+L2-LJ) .lt. 0)        cycle
             if ((L1-L2+LJ) .lt. 0)        cycle
             if ((-L1+L2+LJ) .lt. 0)       cycle
             L=L+1
          end do
       end do
    end do
    allocate(rj(0:LJMAX+1, NRAD), &
         &   ri_mat(NRF*NRF, L))
  end subroutine init_bessel
end module bessel


module xa
  use const, only: R8, C16
  implicit none
  private :: R8, C16
  public; save

  complex(C16), allocatable :: phs(:)
  real(R8),     allocatable :: fj(:,:),dfj(:,:),r(:)
  real(R8)                  :: bk(3),bkrot(3),bkrloc(3)

 contains
   subroutine init_xa(LMAX2,NMAT,NRAD,NB)
     integer, intent(in) :: Lmax2, Nmat, Nrad, NB

     allocate(phs(nb),fj(0:lmax2,nmat),dfj(0:lmax2,nmat))
     allocate(r(nrad))
   end subroutine init_xa
end module xa


module Amn_Mmn
  use const, only: C16

  implicit none
  private; save
  public :: overlap, init_Amn_Mmn, c

  complex(C16), allocatable :: overlap(:,:,:), c(:,:,:)

contains
  subroutine init_Amn_Mmn(Nbands, Npair, Nproj, Nat)
    use param, only: Lmax2

    integer, intent(in) :: Nbands, Npair, Nproj, Nat

    allocate( overlap(Nbands, Nbands, Npair), c(Nproj, (LMAX2+1)**2, Nat) )

    overlap = 0; c = 0
  end subroutine init_Amn_Mmn
end module Amn_Mmn


module pairs
  implicit none
  public; save

  integer, allocatable :: KP(:), KPB(:)
  integer, allocatable :: BQX(:),BQY(:),BQZ(:), BQX1(:),BQY1(:),BQZ1(:)

contains
  subroutine  init_pairs(Npair)
    integer, intent(in) :: Npair

    allocate( kp  (Npair), kpb (Npair), &
         &    bqx (Npair), bqy (Npair), bqz (Npair), &
         &    bqx1(Npair), bqy1(Npair), bqz1(Npair))
  end subroutine init_pairs
end module pairs


module lolog
  use param, only: Lmax2, lomax

  implicit none
  private; save

  integer, public :: Nlo, Nlov, Nlon, n_rad(0:lmax2), ilo(0:lomax)
  logical, public :: loor(0:lomax), lapw(0:lmax2)
end module lolog


module atspdt
  use param, only: Lmax2, Nrf
  use const, only: R8

  implicit none
  private; save

  ! radial function and its slope at RMT
  real(R8), public :: P(0:Lmax2, Nrf), DP(0:Lmax2,Nrf)
end module atspdt


module loabc
  use param, only: lomax, Nloat, Nrf
  use const, only: R8

  implicit none
  private; save

  ! abc calculates the cofficients a,b,c of the lo
  real(R8), public :: alo(0:LOmax, Nloat, Nrf)
end module loabc


module radfu
  use param, only: Nrad, Lmax2, Nrf
  use const, only: R8

  implicit none
  private; save

  ! radial functions large and small component
  real(R8), public :: RF1(Nrad,0:Lmax2,Nrf), RF2(Nrad,0:Lmax2,Nrf)
end module radfu


!---------------------------  Procedure modules  ---------------------------
module     abc_m; contains
subroutine abc(l, rmt, pei, pi12lo, pe12lo, jlo, lapw)
  use param,  only: unit_out, Nrf
  use loabc,  only: alo
  use atspdt, only: P, DP
  use const,  only: R8

  implicit none

  integer,  intent(in) :: l, jlo
  real(R8), intent(in) :: rmt, pei, pi12lo, pe12lo
  logical,  intent(in) :: lapw

  integer  :: irf
  real(R8) :: xac, xbc, clo, alonorm

  character(len=*), parameter :: &
       fmt_alo = "('LO COEFFICIENT: l,A,B,C  ',i2,5X,6F12.5)"

  do irf=1,nrf
     alo(l,jlo,irf)=0.d0
  end do
  if (lapw) then
     irf=2+jlo
     xac=p(l,irf)*dp(l,2)-dp(l,irf)*p(l,2)
     xac= xac * rmt**2
     xbc=p(l,irf)*dp(l,1)-dp(l,irf)*p(l,1)
     xbc= -xbc * rmt**2
     clo=xac*(xac+2.0D0*pi12lo)+xbc* &
          (xbc*pei+2.0D0*pe12lo)+1.0D0
     clo=1.0D0/sqrt(clo)
     write(unit_out,*)clo
     if (clo.gt.2.0D2) clo=2.0d2
     alo(l,jlo,1)=clo*xac
     alo(l,jlo,2)=clo*xbc
     alo(l,jlo,irf)=clo
     write(unit_out, fmt_alo) l, alo(l,jlo,1), alo(l,jlo,2), alo(l,jlo,irf)
  else
     if (jlo.eq.1) then
        alonorm=sqrt(1.d0+(P(l,1)/P(l,2))**2 * PEI)
        alo(l,jlo,1)=1.d0/alonorm
        alo(l,jlo,2)=-P(l,1)/P(l,2)/alonorm
     else
        xbc=-P(l,1)/P(l,1+jlo)
        xac=sqrt(1+xbc**2+2*xbc*PI12LO)
        alo(l,jlo,1)=1.d0/xac
        alo(l,jlo,1+jlo)=xbc/xac
     end if
     write (unit_out, fmt_alo) l, alo(l, jlo, :)
  end if
  return
end subroutine abc
end module     abc_m


module     outwin_m; contains
subroutine outwin(stru, jatom, V, EH, FL, VAL, SLO, Nodes)
  !         Integration der skalarrel. Schroedingergleichung
  !
  !    Rydberg Einheiten

  use param,     only: Nrad
  use const,     only: clight, R8
  use uhelp,     only: A, B
  use structmod, only: struct_t

  implicit none

  !  Input:
  !    stru  struct type
  !    jatom number of the atom
  !    EH    Energie in Hartree
  !    FL    Drehimpuls
  !    V     rad.sym. Potential in Hartree
  type(struct_t), intent(in) :: stru
  real(R8),       intent(in) :: EH, fl, V(Nrad)
  integer,        intent(in) :: jatom

  !  Output:
  !    VAL,SLO:  Wellenfunktion und Steigung am Kugelrand
  !    Nodes:    Anzahl Knoten
  real(R8), intent(out) :: val, slo
  integer,  intent(out) :: nodes

  real(R8) :: D(2,3), Rnet(Nrad), C, E
  real(R8) :: zz, fllp1, s, sf, f0, aa, r, drdi, dg1,dg2,dg3, df1, df2, df3
  real(R8) :: phi, u, x, y, det, b1,b2
  integer  :: i, k

  real(R8), parameter :: H83 = 8/3._R8
  real(R8), parameter :: R83SQ = 64/9._R8
  real(R8), parameter :: R1 = 1/9._R8
  real(R8), parameter :: R2 = -5*R1
  real(R8), parameter :: R3 = 19*R1
  real(R8), parameter :: G0 = 1

  ! Hartree in Ryd
  E = 2*EH

  do i = 1,stru%Npt(jatom)
     Rnet(i) = stru%R0(jatom) * exp(stru%dx(jatom) * (i-1))
  end do

  Nodes = 0
  ZZ = 2*stru%Z(jatom)

  C = merge(2*clight, 1e10_R8, stru%rel)

  FLLP1 = FL*(FL + 1)

  if (stru%Z(jatom) .lt. 0.9D0) then
     S = FL+1.d0
     SF = FL
     F0 = FL/C
  else
     AA = ZZ/C
     S = DSQRT(FLLP1 + 1.D0 - AA*AA)
     SF = S
     F0 = G0*(S - 1.D0)/AA
  endif
  do K = 1,3
     R = RNET(K)
     DRDI = stru%dx(jatom)*R
     A(K) = (R**S)*G0
     B(K) = (R**SF)*F0
     D(1,K) = DRDI*A(K)*S/R
     D(2,K) = DRDI*B(K)*SF/R
  end do

  DG1 = D(1,1)
  DG2 = D(1,2)
  DG3 = D(1,3)
  DF1 = D(2,1)
  DF2 = D(2,2)
  DF3 = D(2,3)
  do K = 4, stru%Npt(jatom)
     R = RNET(K)
     DRDI = stru%dx(jatom)*R

     !       Faktor zwei vor V wegen Hartree-Rydberg !
     PHI = (E - 2.d0*V(K)/R)*DRDI/C
     U = DRDI*C + PHI
     X = -DRDI/R
     Y = -FLLP1*X*X/U + PHI
     DET = R83SQ - X*X + U*Y
     B1 = A(K-1)*H83 + R1*DG1 + R2*DG2 + R3*DG3
     B2 = B(K-1)*H83 + R1*DF1 + R2*DF2 + R3*DF3
     A(K) = (B1*(H83-X) + B2*U)/DET
     B(K) = (B2*(H83+X) - B1*Y)/DET
     if (A(K)*A(K-1) .lt. 0D0) Nodes = Nodes + 1
     DG1 = DG2
     DG2 = DG3
     DG3 = U*B(K) - X*A(K)
     DF1 = DF2
     DF2 = DF3
     DF3 = X*B(K) - Y*A(K)
  end do

  B(1:stru%Npt(jatom))=B(1:stru%Npt(jatom))*c/2

  VAL = A(stru%Npt(jatom))  / RNET(stru%Npt(jatom))
  SLO = DG3/(stru%dx(jatom) * RNET(stru%Npt(jatom)))
  SLO = (SLO-VAL) / RNET(stru%Npt(jatom))
end subroutine outwin
end module outwin_m


module     rint13_m; contains
subroutine rint13(stru, jatom, A, B, X, Y, S)
  !     PERFORM RADIAL INTEGRALS REQUIRED BY BHDK13
  !                            D.D.KOELLING

  use param,     only: Nrad
  use const,     only: clight, R8
  use structmod, only: struct_t

  implicit none

  type(struct_t), intent(in)  :: stru
  integer,        intent(in)  :: jatom
  real(R8),       intent(in)  :: A(Nrad), B(Nrad), X(Nrad), Y(Nrad)
  real(R8),       intent(out) :: S

  integer  :: j, j1
  real(R8) :: d, cin, r,r1, z2,z4, p1,p2

  cin = merge(1/clight**2, 1e-22_R8, stru%rel)

  D=exp(stru%DX(JATOM))

  J=3-mod(stru%Npt(JATOM),2)
  J1=J-1
  R=stru%R0(JATOM)*(D**(J-1))
  R1=R/D
  Z4=0
  Z2=0
10 Z4=Z4+R*(A(J)*X(J)+CIN*B(J)*Y(J))
  R=R*D
  J=J+1
  if(J >= stru%Npt(JATOM)) goto 20
  Z2=Z2+R*(A(J)*X(J)+CIN*B(J)*Y(J))
  R=R*D
  J=J+1
  goto 10
20 P1=stru%R0(JATOM)*(A(1)*X(1)+CIN*B(1)*Y(1))
  P2=R1*(A(J1)*X(J1)+CIN*B(J1)*Y(J1))
  S=2*Z2+4*Z4+R*(A(J)*X(J)+CIN*B(J)*Y(J))+P2
  S=(stru%DX(JATOM)*S+P1)/3.0D0
  if(J1.gt.1) S=S+0.5D0*stru%DX(JATOM)*(P1+P2)
end subroutine rint13
end module     rint13_m


module     atpar_m; contains
subroutine atpar(stru, jatom, itape, jtape)
!!! calculate radial functions for atoms JATOM

  use param,      only: unit_out, Nrad, Nloat, lomax, Lmax2
  use structmod,  only: struct_t
  use lolog,      only: nlo,nlov,nlon,loor,ilo,lapw,n_rad
  use atspdt,     only: P, DP
  use const,      only: R8, clight
  use uhelp,      only: A, B
  use radfu,      only: RF1, RF2

  !! procedure includes
  use diracout_m
  use outwin_m
  use rint13_m
  use dergl_m
  use abc_m

  implicit none

  type(struct_t), intent(in) :: stru
  integer,        intent(in) :: jatom, itape, jtape

  real(R8) :: VR(Nrad), AE(Nrad), BE(Nrad)
  logical  :: rlo(1:nloat, 0:lomax)
  real(r8) :: emist(0:lomax,nloat),E(0:LMAX2),elo(0:LOMAX,nloat),pei(0:lmax2)
  integer  :: imax,irf, jlo, kappa, i,k,l, node,nodes,nodel, m
  real(R8) :: dele,delei, fl, ei,e1, uvb,duvb,uv,duv,uve,duve, ovlp, trx
  real(R8) :: try, r_m, pi12lo, pe12lo, cross

  !.....READ TOTAL SPHERICAL POTENTIAL V(0,0) OF TAPEjtape=VSP
  !     NORM OF V(0,0)=V00(R)*R/SQRT(4.D0*PI)
  read(jtape, '(3X)')
  read(jtape, '(15X,I3//)') i ! i is not used, but apparently necessary for
  read(jtape, '(/)')          ! correct reading
  read(jtape, '(3X,4E19.12)') VR(1 : stru%Npt(jatom))
  read(jtape, '(/)')
  read(jtape, '(///)')

  VR(1:stru%Npt(jatom)) = VR(1:stru%Npt(jatom)) / 2

  write(unit_out,*)'ATPAR'
  nlo=0
  nlov=0
  nlon=0
  ilo=0
  n_rad=2

  atoms: do I=1,JATOM
     read(itape)E
     read(itape)elo
     if(i == jatom) then
        do l=0,lmax2
           lapw(l)=.true.
           if(e(l) > 150) then
              e(l)=e(l)-200
              lapw(l)=.false.
           endif
        enddo
     endif
     LOs: do l = 0,lomax
        loor(l)=.false.
        do k=1,nloat
           rlo(k,l)=.false.
           if (i == jatom) then
              if (elo(l,k) < 995) then
                 ilo(l)=ilo(l)+1
                 nlo = nlo + (2*l+1) * stru%mult(i)
                 if(.not. lapw(l).and.k == 1) cycle
                 if(k == nloat) then
                    rlo(ilo(l),l)=.true.
                    cycle
                 endif
                 loor(l)=.true.
              endif
           else
              if (elo(l,k) < 995) nlov = nlov + (2*l+1)*stru%mult(i)
           endif
        end do
     end do LOs
  end do atoms

  if(jatom /= stru%Nneq) then
     do I=JATOM+1,stru%Nneq
        read(itape) EMIST
        read(itape) EMIST
        do l=0,lomax
           do k=1,nloat
              if (emist(l,k) < 995)  &
                   nlon = nlon + (2*l+1)*stru%mult(i)
           end do
        end do
     end do
  end if

  write(unit_out, "(/10X,'ATOMIC PARAMETERS FOR ',A10/)") stru%aname(JATOM)
  write(unit_out, "(10X,' ENERGY PARAMETERS ARE',7F7.2)") E
  write(unit_out, "(/11X,1HL,5X,4HU(R),10X, 5HU'(R),9X,5HDU/DE,8X,6HDU'/DE,6X,7HNORM-U')")

  lloop: do l=0,LMAX2
     DELE=2.0D-3
     DELEI=0.25D0/DELE
     FL=L
     EI=E(l)/2.0d0
     !     CALCULATE ENERGY-DERIVATIVE BY FINITE DIFFERENCE
     !     DELE IS THE UPWARD AND DOWNWARD ENERGY SHIFT IN HARTREES

     E1=EI-DELE
     ! outwin() sets A(:), B(:)!
     call outwin(stru, jatom, Vr, E1, FL, UVB, DUVB, NODEL)
     call rint13(stru, jatom, A, B, A, B, OVLP)
     TRX=1.0D0/sqrt(OVLP)
     IMAX=stru%Npt(JATOM)
     do M=1,IMAX
        AE(M)=TRX*A(M)
        BE(M)=TRX*B(M)
     end do
     UVB=TRX*UVB
     DUVB=TRX*DUVB
     E1=EI+DELE
     call outwin(stru, jatom, Vr, E1, FL, UVE, DUVE, NODE)
     call rint13(stru, jatom, A, B, A, B, OVLP)
     TRX=1.0d0/sqrt(OVLP)
     UVE=DELEI*(TRX*UVE-UVB)
     DUVE=DELEI*(TRX*DUVE-DUVB)
     IMAX=stru%Npt(JATOM)
     do M=1,IMAX
        AE(M)=DELEI*(TRX*A(M)-AE(M))
        BE(M)=DELEI*(TRX*B(M)-BE(M))
     end do

     !     Calculate function at EI
     call outwin(stru, jatom, Vr, EI, FL, UV, DUV, NODES)
     call rint13(stru, jatom, A, B, A, B, OVLP)
     TRX=1.0d0/sqrt(OVLP)
     P(l,1)=TRX*UV
     DP(l,1)=TRX*DUV
     IMAX=stru%Npt(JATOM)
     do M=1,IMAX
        A(M)=TRX*A(M)
        B(M)=TRX*B(M)
     end do

     !     Ensure orthogonalization
     call rint13(stru, jatom, A, B, AE, BE, CROSS)
     TRY=-CROSS
     IMAX=stru%Npt(JATOM)
     do M=1,IMAX
        AE(M)=(AE(M)+TRY*A(M))
        BE(M)=(BE(M)+TRY*B(M))
     end do
     IMAX=stru%Npt(JATOM)
     do I=1,IMAX
        RF1(I,l,1)=A(I)
        RF2(I,l,1)=B(I)
        RF1(I,l,2)=AE(I)
        RF2(I,l,2)=BE(I)
     end do
     P(l,2)=UVE+TRY*P(l,1)
     DP(l,2)=DUVE+TRY*DP(l,1)
     call RINT13(stru, jatom, AE, BE, AE, BE, PEI(l))
     write(unit_out, "(10X,I2,5E14.6,5X,3I2)") &
          L,P(l,1),DP(l,1),P(l,2),DP(l,2)
  end do lloop
!
! nun fur lo
!
  loloop: do l=0,lomax
     irf=2
     iloloop: do jlo=1,ilo(l)
        if (lapw(l) .or. jlo/=1) then
           irf=irf+1
           DELE=2.0D-3
           DELEI=0.25D0/DELE
           FL=L
           EI=elo(l,jlo)/2.d0
           !
           !     CALCULATE FUNCTION AT EI
           if(rlo(jlo,l)) then
              ei=elo(l,nloat)/2.d0
              kappa=l
              call diracout(stru, jatom, Vr, ei, kappa, uv, duv, nodes)
              call dergl(stru, jatom, a, b)
              do m = 1, stru%Npt(jatom)
                 r_m = stru%R0(jatom) * exp(stru%dx(jatom) * (m-1))
                 b(m) = b(m)*r_m/(2.d0*clight+(elo(l,jlo)- &
                      2.d0*vr(m)/r_m)/(2.d0*clight))
                 b(m)=b(m)*clight
              enddo
           else
              call outwin(stru, jatom, Vr, ei, fl, uv, duv, nodes)
           endif

           call RINT13(stru, jatom, A, B, A, B, OVLP)
           TRX=1.0d0/sqrt(OVLP)
           P(l,irf)=TRX*UV
           DP(l,irf)=TRX*DUV
           IMAX=stru%Npt(JATOM)
           n_rad(l)=irf
           do M=1,IMAX
              rf1(M,l,irf)=TRX*A(M)
              rf2(M,l,irf)=TRX*B(M)
           end do

           call RINT13(stru, jatom, rf1(1,l,1), rf2(1,l,1), &
                &      rf1(1,l,irf), rf2(1,l,irf), pi12lo)
           call RINT13(stru, jatom, rf1(1,l,2), rf2(1,l,2), &
                &      rf1(1,l,irf), rf2(1,l,irf), pe12lo)
        end if
        call abc (l, stru%RMT(jatom), pei(l), pi12lo, pe12lo, jlo, lapw(l))
     end do iloloop
  end do loloop

  write(unit_out, "('number of rad. functions per L:',8I3)") &
       n_rad(0 : lmax2)
end subroutine atpar
end module     atpar_m


module                 gaunt1_m; contains
real(R8) pure function gaunt1(LP,L,LS,MP,M,MS)
!
  use assleg, only: YR, N
  use const,  only: R8

  implicit none
!
!        Arguments
!
  integer, intent(in) :: L, LP, LS, M, MP, MS
!
!     ..................................................................
!
!        GAUNT computes the integral of
!           CONJG(Y(LP,MP))*Y(L,M)*Y(LS,MS) for LP+L+LS  <=  23
!        using gaussian quadrature with N=12 as given by
!           M. Abramowitz and I.A. Stegun,
!           'Handbook of Mathematical Functions',
!           NBS Applied Mathematics Series 55 (1968), pages 887 and 916
!        written by Bruce Harmon based on suggestion by W. Rudge
!        Iowa State Sept.1973
!
!        extended by M. Weinert and E. Wimmer
!        Northwestern University March 1980
!
!     ..................................................................
!
!      INTEGER            MAXDIM, N
!      PARAMETER          (MAXDIM = 81, N = 6)
!
!        Local Scalars
!
      integer  :: I, IL, ILP, ILS
      real(R8) :: S
      real(R8), parameter :: W(N) = (/ &
           0.24914704581340D+0, 0.23349253653836D+0, &
           0.20316742672307D+0, 0.16007832854335D+0, &
           0.10693932599532D+0, 0.04717533638651D+0 /)

      IL = L*(L+1) + M + 1
      ILP = LP*(LP+1) + MP + 1
      ILS = LS*(LS+1) + MS + 1
      S = 0
      do I = 1, N
         S = S + W(I)*YR(I,ILP)*YR(I,IL)*YR(I,ILS)
      end do
      GAUNT1 = S
!
end function gaunt1
end module   gaunt1_m


module     gaunt2_m; contains
subroutine gaunt2
  ! set YR needed in function GAUNT1

  use assleg, only: YR, N, init_assleg
  use const,  only: R8
  implicit none

  real(R8), parameter :: SNULL = 1.0e-10_R8
  real(R8), parameter :: X(N) = &
       (/ 0.12523340851147D+0, 0.36783149899818D+0, &
       &  0.58731795428662D+0, 0.76990267419431D+0, &
       &  0.90411725637048D+0, 0.98156063424672D+0  /)

  integer  :: I, IDWN, IX, K, L, L1, L2, LM, LM2, LOMAX, M, M1
  real(R8) :: C1L, C2L, CD, CSR, CTH, CYP, FACTOR, FPI, RF, SGNM, STH, TCTH
  real(R8) :: P(10,10)

  call init_assleg

  FPI = 16*atan(1._R8)
  FACTOR = FPI**(1/3._R8)
  LOMAX = 8
!
  do K = 1, N
     CTH = X(K)
     STH = sqrt(1-CTH*CTH)
     RF = 1/sqrt(FPI)
     YR(K,1) = RF*FACTOR
     I = 1
     P(1,1) = 1
     P(2,1) = CTH
     C2L = CTH
     TCTH = CTH + CTH
     L1 = 2
     L = 1
10   continue
     M = 1
     I = I + L
     IDWN = I + 2
     M1 = 2
     L2 = L
     L = L1
     L1 = L + 1
     LM = L2
     LM2 = L
     CD = 1
     C2L = C2L + TCTH
     SGNM = 1
20   continue
     !
     !        recurse upward in L
     !
     P(L1,M) = (C2L*P(L,M)-LM*P(L2,M))/LM2
     C1L = (LM+1)*CTH
     P(L,M1) = 0
     !
     !        recurse upward in M
     !
     if (abs(STH)  >=  SNULL) P(L,M1) = (C1L*P(L,M)-LM2*P(L1,M))/STH
30   continue
     I = I + 1
     IDWN = IDWN - 1
     CSR = sqrt((2*L-1)/(FPI*CD))
     CYP = SGNM*CSR*P(L,M)
     YR(K,I) = CYP*FACTOR
     IX = I - (L*L - L + 1)
     if (IDWN .ne. I) YR(K,IDWN) = FACTOR*CYP*(-1)**IX
     M = M1
     M1 = M + 1
     LM2 = LM2 - 1
     LM = LM + 1
     CD = CD*LM*LM2
     SGNM = -SGNM
     if (M - L < 0) goto 20
     if (M - L ==0) goto 30
     if (L <= LOMAX) goto 10
  end do
end subroutine gaunt2
end module gaunt2_m


module     harmon_m; contains
subroutine harmon(stru, jatom, N, X, Y, Z, Lmax2, F, DF)
  use const,     only: R8
  use structmod, only: struct_t

  !! procedure includes
  use dvbes1_m
  use sphbes_m

  implicit none

  type(struct_t), intent(in)  :: stru
  integer,        intent(in)  :: N, Lmax2, jatom
  real(R8),       intent(in)  :: X(N), Y(N), Z(N)
  real(R8),       intent(out) :: F(Lmax2+1,N), DF(Lmax2+1,N)

  real(R8) :: A(3), xm, xa
  integer  :: LMX, i, j

  LMX=Lmax2+1
  do I=1,N
     A = matmul(stru%gbas, (/X(I), Y(I), Z(I)/))

     XM = sqrt(sum(A**2))
     XA = stru%RMT(jatom) * XM
     call SPHBES(Lmax2, XA, F(:,I))
     call DVBES1(F(:,I), DF(:,I), XA, LMX)
     do J=1,LMX
        DF(J,I)=XM*DF(J,I)
     end do
  end do
end subroutine harmon
end module harmon_m


module     radint_m; contains
subroutine radint(stru, jatom, ljmax, bm)
  use param,     only: Lmax2, Nrad
  use bessel,    only: rj, ri_mat
  use lolog,     only: n_rad
  use const,     only: R8
  use radfu,     only: RF1, RF2
  use structmod, only: struct_t

  !! procedure includes
  use sphbes_m
  use rint13_m

  implicit none

  type(struct_t), intent(in) :: stru
  integer,        intent(in) :: jatom, LJmax
  real(R8),       intent(in) :: bm

  real(R8) :: A(Nrad), B(Nrad), X(Nrad), Y(Nrad), RX
  integer  :: L_index,l1,l2,lj, R_index, i, if1,if2

  do  I=1,stru%Npt(JATOM)
     RX=stru%R0(JATOM)*exp(stru%dx(JATOM)*(i-1))*BM
     call sphbes(LJMAX+1,RX,rj(:,I))
  enddo
  L_index=0
  do L1=0,LMAX2
     do L2=0,LMAX2
        ljloop: do LJ=0,LJMAX
           if (mod((L1+L2+LJ),2) == 1  &
                .or. (L1+L2-LJ)   < 0  &
                .or. (L1-L2+LJ)   < 0  &
                .or. (-L1+L2+LJ)  < 0) &
                cycle ljloop
           L_index=L_index+1
           R_index=0
           do IF1=1,n_rad(l1)
              do IF2=1,n_rad(l2)
                 R_index=R_index+1
                 do  I=1,stru%Npt(JATOM)
                    A(i)=rf1(i,l1,if1)*rj(lj,i)
                    B(i)=rf2(i,l1,if1)*rj(lj,i)
                    X(i)=rf1(i,l2,if2)
                    Y(i)=rf2(i,l2,if2)
                 enddo
                 call RINT13(stru, jatom, A,B, X,Y, ri_mat(r_index,l_index))

              end do
           end do
        end do ljloop
     end do
  end do
end subroutine radint
end module radint_m


!!/---
!! Local Variables:
!! mode: f90
!! End:
!!\---
!!
!! Time-stamp: <2016-08-02 13:22:26 assman@faepop71.tu-graz.ac.at>
