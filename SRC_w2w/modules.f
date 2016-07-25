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
!! gener:   BR1(3,3), BR2(3,3)
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
!! struct:  init_struct(), nat, iord, lattic, cform, title, irel, rel,
!!          AA,BB,CC, VOL, pia(3), alpha(3), aname(:), R0(:), DX(:),
!!          RMT(:), zz(:), mult(:), jrj(:), v(:), pos(:,:), trans(:,:),
!!          transij(:,:), rotloc(:,:,:), rotij(:,:,:), iz(:,:,:)
!!
!! w2w:     unit_amn=7,  unit_mmn  =8, unit_nnkp=11, unit_eig=12,
!!          unit_ene=50, unit_fermi=51, iBlock=128, NMAT, NDIF
!!
!! xa:      init_xa(), bk(3), bkrot(3), bkrloc(3), fj(:,:), dfj(:,:),
!!          phs(:), r(:)
!!
!! Procedure modules:
!!
!!    abc_m, atpar_m, gaunt1_m, gaunt2_m, harmon_m, latgen_m,
!!    outwin_m, radint_m, rint13_m, rotate_m, rotdef_m
!!
!!\===============================================


module w2w
  implicit none
  public

  integer, parameter :: unit_amn=7,  unit_mmn  =8, unit_nnkp=11, unit_eig=12
  integer, parameter :: unit_ene=50, unit_fermi=51
  ! Optimize IBLOCK for your hardware (32-255)
  integer, parameter :: iBlock=128

  integer :: NMAT=0, NDIF=0
end module w2w


module assleg
  use const, only: R8

  implicit none
  private

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


module struct
  use const, only: R8
  implicit none
  private :: R8

  public

  logical                    :: rel
  real(R8)                   :: AA,BB,CC,VOL,pia(3),alpha(3)
  real(R8),allocatable       :: R0(:),DX(:),RMT(:),zz(:),rotloc(:,:,:),v(:)
  real(R8),allocatable       :: trans(:,:)
  real(R8),allocatable       :: pos(:,:)
  real(R8),allocatable       :: rotij(:,:,:),transij(:,:)
  character(4)               :: lattic,irel,cform
  character(80)              :: title
  character(10), allocatable :: aname(:)
  integer                    :: nat,iord
  integer,allocatable        :: mult(:),jrj(:)
  integer,allocatable        :: iz(:,:,:)

 contains
  subroutine init_struct
    use param,      only: unit_struct, unit_out
    use w2w,        only: Ndif
    use reallocate, only: realloc
    use const,      only: R8

    implicit none

    integer  :: ios, inum, isplit, iatnr
    real(R8) :: test
    integer  :: index, i, j, j1, j2, m, jatom

    parameter (test=1.e-5_R8)

    read (unit_struct,1000) title
    read (unit_struct,1010) lattic,nat,cform,irel
    REL=.true.
    if(IREL.eq.'NREL') REL=.false.
    allocate(aname(nat),mult(nat),jrj(nat),r0(nat),dx(nat),rmt(nat))
    allocate(zz(nat),rotloc(3,3,nat),v(nat))
    v=0.0d0
    allocate (pos(3,48*nat))
    read (unit_struct,1020) aa,bb,cc,alpha(1),alpha(2),alpha(3)
    if(abs(ALPHA(1)).lt.test) ALPHA(1)=90
    if(abs(ALPHA(2)).lt.test) ALPHA(2)=90
    if(abs(ALPHA(3)).lt.test) ALPHA(3)=90
    INDEX=0
    do jatom=1,NAT
     INDEX=INDEX+1
     read(unit_struct,1030,iostat=ios) &
          iatnr, (pos(j,index),j=1,3), mult(jatom), isplit
     if(ios /= 0) then
      write(unit_out,*) iatnr,(pos(j,index),j=1,3),mult(jatom)
      write(unit_out,*) 'ERROR IN STRUCT FILE READ'
      stop
      endif
      if (mult(jatom) .eq. 0) then
       write (unit_out,6000) jatom, index, mult(jatom)
       stop
       endif
       do m=1,mult(jatom)-1
        index=index+1
        read(unit_struct,1031)iatnr,(pos(j,index),j=1,3)
       enddo
       read(unit_struct,1050) aname(jatom),jrj(jatom),r0(jatom),rmt(jatom), &
            zz(jatom)
       dx(jatom)=log(rmt(jatom)/r0(jatom)) / (jrj(jatom)-1)
       rmt(jatom)=r0(jatom)*exp( dx(jatom)*(jrj(jatom)-1) )
       read(unit_struct,1051) ((rotloc(i,j,jatom),i=1,3),j=1,3)
    enddo
    ndif=index
    call realloc(pos, (/ 3, ndif /))
    allocate(rotij(3,3,ndif),transij(3,ndif))
    read(unit_struct,1151)iord
    allocate(iz(3,3,iord),trans(3,iord))
    do j=1,iord
       read(unit_struct,1101)((iz(j1,j2,j),j1=1,3),trans(j2,j),j2=1,3),inum
    enddo

1000 format(A80)
1010 format(A4,23X,I3,1x,a4,/,13X,A4,18X,A4)
1020 format(6F10.7,10X,F10.7)
1030 format(4X,I4,4X,F10.7,3X,F10.7,3X,F10.7,/,15X,I2,17X,I2)
1031 format(4X,I4,4X,F10.7,3X,F10.7,3X,F10.7)
1050 format(A10,5X,I5,5X,F10.9,5X,F10.5,5X,F10.5)
1051 format(20X,3F10.8)
1101 format(3(3I2,F11.8/),I8)
1151 format(I4)
6000 format(///,3X,'ERROR IN LAPW0 : MULT(JATOM)=0 ...', &
          /, 20X,'JATOM=',I3,3X,'INDEX=',I3,3X,'MULT=',I3)
  end subroutine init_struct
end module struct


module bessel
  use const, only: R8
  implicit none
  private :: R8

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
  public

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
  private
  public :: overlap, init_Amn_Mmn, c

  complex(C16), allocatable :: overlap(:,:,:), c(:,:,:)

contains
  subroutine init_Amn_Mmn(nbands,npair,nproj)
    use param, only: Lmax2
    use w2w,   only: Ndif

    integer, intent(in) :: nbands,npair,nproj

    allocate(overlap(nbands,nbands,npair), c(nproj,(LMAX2+1)*(LMAX2+1),ndif))

    overlap = 0; c = 0
  end subroutine init_Amn_Mmn
end module Amn_Mmn


module pairs
  implicit none
  public

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
  private

  integer, public :: Nlo, Nlov, Nlon, n_rad(0:lmax2), ilo(0:lomax)
  logical, public :: loor(0:lomax), lapw(0:lmax2)
end module lolog


module gener
  use const, only: R8

  implicit none
  private

  ! transformation between u.c. and cartesian coordinates
  real(R8), public :: BR1(3,3), BR2(3,3)
end module gener


module atspdt
  use param, only: Lmax2, Nrf
  use const, only: R8

  implicit none
  private

  ! radial function and its slope at RMT
  real(R8), public :: P(0:Lmax2, Nrf), DP(0:Lmax2,Nrf)
end module atspdt


module loabc
  use param, only: lomax, Nloat, Nrf
  use const, only: R8

  implicit none
  private

  ! abc calculates the cofficients a,b,c of the lo
  real(R8), public :: alo(0:LOmax, Nloat, Nrf)
end module loabc


module radfu
  use param, only: Nrad, Lmax2, Nrf
  use const, only: R8

  implicit none
  private

  ! radial functions large and small component
  real(R8), public :: RF1(Nrad,0:Lmax2,Nrf), RF2(Nrad,0:Lmax2,Nrf)
end module radfu


!---------------------------  Procedure modules  ---------------------------
module     abc_m; contains
subroutine abc(l,jatom,pei,pi12lo,pe12lo,jlo,lapw)
  use param,  only: unit_out, Nrf
  use struct, only: RMT
  use loabc,  only: alo
  use atspdt, only: P, DP
  use const,  only: R8

  implicit none

  integer,  intent(in) :: l, jatom, jlo
  real(R8), intent(in) :: pei, pi12lo, pe12lo
  logical,  intent(in) :: lapw

  integer  :: irf, jrf
  real(R8) :: xac, xbc, clo, alonorm
  !---------------------------------------------------------------------
  !
  do irf=1,nrf
     alo(l,jlo,irf)=0.d0
  end do
  if (lapw) then
     irf=2+jlo
     xac=p(l,irf)*dp(l,2)-dp(l,irf)*p(l,2)
     xac=xac*rmt(jatom)*rmt(jatom)
     xbc=p(l,irf)*dp(l,1)-dp(l,irf)*p(l,1)
     xbc=-xbc*rmt(jatom)*rmt(jatom)
     clo=xac*(xac+2.0D0*pi12lo)+xbc* &
          (xbc*pei+2.0D0*pe12lo)+1.0D0
     clo=1.0D0/sqrt(clo)
     write(unit_out,*)clo
     if (clo.gt.2.0D2) clo=2.0d2
     alo(l,jlo,1)=clo*xac
     alo(l,jlo,2)=clo*xbc
     alo(l,jlo,irf)=clo
     write(unit_out,10) l, alo(l,jlo,1), alo(l,jlo,2), alo(l,jlo,irf)
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
     write (unit_out,10)l,(alo(l,jlo,jrf),jrf=1,nrf)
  end if
10 format ('LO COEFFICIENT: l,A,B,C  ',i2,5X,6F12.5)
  return
end subroutine abc
end module     abc_m


module     outwin_m; contains
subroutine outwin(REL,V,RNOT,DH,JRI,EH,FL,VAL,SLO,Nodes,Z)
  !         Integration der skalarrel. Schroedingergleichung
  !
  !    Rydberg Einheiten

  use param, only: Nrad
  use const, only: clight, R8
  use uhelp, only: A, B

  implicit none

  real(R8) :: EH, fl, Z, V(Nrad), Rnot, dh, val, slo
  integer  :: jri, nodes
  logical  :: rel

  !  Input:
  !    EH    Energie in Hartree
  !    FL    Drehimpuls
  !    Z     Kernladung
  !    V     rad.sym. Potential in Hartree
  !    RNOT  erster radialer Netzpunkt
  !    DH    log. Schrittweite
  !    JRI   Anzahl radialer Netzpunkte
  intent(in)  :: EH, fl, Z, V, Rnot, dh, jri

  !  Output:
  !    VAL,SLO:  Wellenfunktion und Steigung am Kugelrand
  !    Nodes:    Anzahl Knoten
  intent(out) :: val, slo, nodes

  real(R8) :: D(2,3), Rnet(Nrad), C, E
  real(R8) :: zz, fllp1, s, sf, f0, aa, r, drdi, dg1,dg2,dg3, df1, df2, df3
  real(R8) :: phi, u, x, y, det, b1,b2
  integer  :: iiij, k

  real(R8), parameter :: H83 = 8/3._R8
  real(R8), parameter :: R83SQ = 64/9._R8
  real(R8), parameter :: R1 = 1/9._R8
  real(R8), parameter :: R2 = -5*R1
  real(R8), parameter :: R3 = 19*R1
  real(R8), parameter :: G0 = 1

  ! Hartree in Ryd
  E = 2*EH

  do iiij=1,JRI
     RNET(iiij)=RNOT*(exp(DH*(iiij-1)))
  end do

  Nodes = 0
  ZZ = Z + Z

  C = merge(2*clight, 1e10_R8, rel)

  FLLP1 = FL*(FL + 1)

  IF (Z .LT. 0.9D0) THEN
     S = FL+1.d0
     SF = FL
     F0 = FL/C
  ELSE
     AA = ZZ/C
     S = DSQRT(FLLP1 + 1.D0 - AA*AA)
     SF = S
     F0 = G0*(S - 1.D0)/AA
  ENDIF
  do K = 1,3
     R = RNET(K)
     DRDI = DH*R
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
  do K = 4, JRI
     R = RNET(K)
     DRDI = DH*R

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
     IF (A(K)*A(K-1) .LT. 0D0) Nodes = Nodes + 1
     DG1 = DG2
     DG2 = DG3
     DG3 = U*B(K) - X*A(K)
     DF1 = DF2
     DF2 = DF3
     DF3 = X*B(K) - Y*A(K)
  end do

  do iiij=1,JRI
     B(iiij)=B(iiij)*c/2.d0
  end do
  !
  VAL = A(JRI)/RNET(JRI)
  SLO = DG3/(DH*RNET(JRI))
  SLO = (SLO-VAL)/RNET(JRI)
end subroutine outwin
end module outwin_m


module     rint13_m; contains
subroutine rint13(A,B,X,Y,S,JATOM)
  !     PERFORM RADIAL INTEGRALS REQUIRED BY BHDK13
  !                            D.D.KOELLING

  use param,  only: Nrad
  use struct, only: rel, dx, jrj, R0
  use const,  only: clight, R8

  implicit none

  integer,  intent(in)  :: jatom
  real(R8), intent(in)  :: A(Nrad), B(Nrad), X(Nrad), Y(Nrad)
  real(R8), intent(out) :: S

  integer  :: j, j1
  real(R8) :: d, cin, r,r1, z2,z4, p1,p2

  cin = merge(1/clight**2, 1e-22_R8, rel)

  D=EXP(DX(JATOM))

  J=3-MOD(JRJ(JATOM),2)
  J1=J-1
  R=R0(JATOM)*(D**(J-1))
  R1=R/D
  Z4=0
  Z2=0
10 Z4=Z4+R*(A(J)*X(J)+CIN*B(J)*Y(J))
  R=R*D
  J=J+1
  IF(J.GE.JRJ(JATOM)) GOTO 20
  Z2=Z2+R*(A(J)*X(J)+CIN*B(J)*Y(J))
  R=R*D
  J=J+1
  GOTO 10
20 P1=R0(JATOM)*(A(1)*X(1)+CIN*B(1)*Y(1))
  P2=R1*(A(J1)*X(J1)+CIN*B(J1)*Y(J1))
  S=2*Z2+4*Z4+R*(A(J)*X(J)+CIN*B(J)*Y(J))+P2
  S=(DX(JATOM)*S+P1)/3.0D0
  IF(J1.GT.1) S=S+0.5D0*DX(JATOM)*(P1+P2)
end subroutine rint13
end module     rint13_m


module     atpar_m; contains
subroutine atpar(JATOM, itape, jtape)
!!! calculate radial functions for atoms JATOM

  use param,  only: unit_out, Nrad, Nloat, lomax, Lmax2
  use struct, only: JRJ, mult, Nat, aname, R0, dx, zz, rel
  use lolog,  only: nlo,nlov,nlon,loor,ilo,lapw,n_rad
  use atspdt, only: P, DP
  use const,  only: R8, clight
  use uhelp,  only: A, B
  use radfu,  only: RF1, RF2

  !! procedure includes
  use diracout_m
  use outwin_m
  use rint13_m
  use dergl_m
  use abc_m

  implicit none

  integer, intent(in) :: jatom, itape, jtape

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
  read(jtape, '(3X,4E19.12)') VR(1 : JRJ(jatom))
  read(jtape, '(/)')
  read(jtape, '(///)')

  VR(1:JRJ(jatom)) = VR(1:JRJ(jatom)) / 2

  write(unit_out,*)'ATPAR'
  nlo=0
  nlov=0
  nlon=0
  ilo=0
  n_rad=2

  atoms: do I=1,JATOM
     read(itape)E
     read(itape)elo
     if(i.eq.jatom) then
        do l=0,lmax2
           lapw(l)=.true.
           if(e(l).gt.150.) then
              e(l)=e(l)-200.d+0
              lapw(l)=.false.
           endif
        enddo
     endif
     LOs: do l = 0,lomax
        loor(l)=.false.
        do k=1,nloat
           rlo(k,l)=.false.
           if (i.eq.jatom) then
              if (elo(l,k).lt.(995.d+0)) then
                 ilo(l)=ilo(l)+1
                 nlo=nlo+((2*l+1))*mult(i)
                 if(.not.lapw(l).and.k.eq.1) cycle
                 if(k.eq.nloat) then
                    rlo(ilo(l),l)=.true.
                    cycle
                 endif
                 loor(l)=.true.
              endif
           else
              if (elo(l,k).lt.(995.d+0)) nlov=nlov+((2*l+1))*mult(i)
           endif
        end do
     end do LOs
  end do atoms

  if(jatom /= nat) then
     do I=JATOM+1,NAT
        read(itape) EMIST
        read(itape) EMIST
        do l=0,lomax
           do k=1,nloat
              if (emist(l,k).lt.(995.0D+0))  &
                   nlon=nlon+((2*l+1))*mult(i)
           end do
        end do
     end do
  end if

  write(unit_out, "(/10X,'ATOMIC PARAMETERS FOR ',A10/)") ANAME(JATOM)
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
     call OUTWIN(REL,VR,R0(JATOM),DX(JATOM),JRJ(JATOM),E1,            &
          FL,UVB,DUVB,NODEL,ZZ(jatom))
     call RINT13(A,B,A,B,OVLP,JATOM)
     TRX=1.0D0/sqrt(OVLP)
     IMAX=JRJ(JATOM)
     do M=1,IMAX
        AE(M)=TRX*A(M)
        BE(M)=TRX*B(M)
     end do
     UVB=TRX*UVB
     DUVB=TRX*DUVB
     E1=EI+DELE
     call OUTWIN(REL,VR,R0(JATOM),DX(JATOM),JRJ(JATOM),E1,            &
          FL,UVE,DUVE,NODE,ZZ(jatom))
     call RINT13(A,B,A,B,OVLP,JATOM)
     TRX=1.0d0/sqrt(OVLP)
     UVE=DELEI*(TRX*UVE-UVB)
     DUVE=DELEI*(TRX*DUVE-DUVB)
     IMAX=JRJ(JATOM)
     do M=1,IMAX
        AE(M)=DELEI*(TRX*A(M)-AE(M))
        BE(M)=DELEI*(TRX*B(M)-BE(M))
     end do
     !
     !     CALCULATE FUNCTION AT EI
     !
     call OUTWIN(REL,VR(1),R0(JATOM),DX(JATOM),JRJ(JATOM),EI,         &
          FL,UV,DUV,NODES,ZZ(jatom))
     call RINT13(A,B,A,B,OVLP,JATOM)
     TRX=1.0d0/sqrt(OVLP)
     P(l,1)=TRX*UV
     DP(l,1)=TRX*DUV
     IMAX=JRJ(JATOM)
     do M=1,IMAX
        A(M)=TRX*A(M)
        B(M)=TRX*B(M)
     end do
     !
     !     INSURE ORTHOGONALIZATION
     !
     call RINT13(A,B,AE,BE,CROSS,JATOM)
     TRY=-CROSS
     IMAX=JRJ(JATOM)
     do M=1,IMAX
        AE(M)=(AE(M)+TRY*A(M))
        BE(M)=(BE(M)+TRY*B(M))
     end do
     IMAX=JRJ(JATOM)
     do I=1,IMAX
        RF1(I,l,1)=A(I)
        RF2(I,l,1)=B(I)
        RF1(I,l,2)=AE(I)
        RF2(I,l,2)=BE(I)
     end do
     P(l,2)=UVE+TRY*P(l,1)
     DP(l,2)=DUVE+TRY*DP(l,1)
     call RINT13(AE,BE,AE,BE,PEI(l),JATOM)
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
              call diracout(rel,vr(1),r0(jatom),dx(jatom),jrj(jatom),    &
                   &        ei,kappa,uv,duv,nodes,zz(jatom))
              call dergl(a,b,r0(jatom),dx(jatom),jrj(jatom))
              do m = 1, jrj(jatom)
                 r_m = r0(jatom)*exp(dx(jatom)*(m-1))
                 b(m) = b(m)*r_m/(2.d0*clight+(elo(l,jlo)- &
                      2.d0*vr(m)/r_m)/(2.d0*clight))
                 b(m)=b(m)*clight
              enddo
           else
              call outwin(rel,vr(1),r0(jatom),dx(jatom),jrj(jatom),   &
                   ei,fl,uv,duv,nodes,zz(jatom))
           endif

           call RINT13(A,B,A,B,OVLP,JATOM)
           TRX=1.0d0/sqrt(OVLP)
           P(l,irf)=TRX*UV
           DP(l,irf)=TRX*DUV
           IMAX=JRJ(JATOM)
           n_rad(l)=irf
           do M=1,IMAX
              rf1(M,l,irf)=TRX*A(M)
              rf2(M,l,irf)=TRX*B(M)
           end do

           call RINT13(rf1(1,l,1),rf2(1,l,1), &
                rf1(1,l,irf),rf2(1,l,irf),pi12lo,JATOM)
           call RINT13(rf1(1,l,2),rf2(1,l,2), &
                rf1(1,l,irf),rf2(1,l,irf),pe12lo,JATOM)
        end if
        call abc (l,jatom,pei(l),pi12lo,pe12lo,jlo,lapw(l))
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
!           CONJG(Y(LP,MP))*Y(L,M)*Y(LS,MS) for LP+L+LS .LE. 23
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
     if (abs(STH) .ge. SNULL) P(L,M1) = (C1L*P(L,M)-LM2*P(L1,M))/STH
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
     if (L .le. LOMAX) goto 10
  end do
end subroutine gaunt2
end module gaunt2_m


module     harmon_m; contains
subroutine harmon(N,X,Y,Z,LMAX2,F,DF,RI)
  use const, only: R8
  use gener, only: br1

  !! procedure includes
  use dvbes1_m
  use sphbes_m

  implicit none

  integer,  intent(in)  :: N, Lmax2
  real(R8), intent(in)  :: X(N),Y(N),Z(N), RI
  real(R8), intent(out) :: F(LMAX2+1,N), DF(LMAX2+1,N)

  real(R8) :: A(3), xm, xa
  integer  :: LMX, i, j

  LMX=LMAX2+1
  do I=1,N
     A(:)=X(I)*BR1(1,:)+Y(I)*BR1(1,:)+Z(I)*BR1(1,:)

     XM=sqrt(A(1)**2 + A(2)**2 + A(3)**2)
     XA=RI*XM
     call SPHBES(LMAX2, XA, F(:,I))
     call DVBES1(F(:,I), DF(:,I), XA, LMX)
     do J=1,LMX
        DF(J,I)=XM*DF(J,I)
     end do
  end do
end subroutine harmon
end module harmon_m


module     radint_m; contains
subroutine radint(JATOM,LJMAX,BM)
  use param,  only: Lmax2, Nrad
  use struct, only: jrj, R0, dx
  use bessel, only: rj, ri_mat
  use lolog,  only: n_rad
  use const,  only: R8
  use radfu,  only: RF1, RF2

  !! procedure includes
  use sphbes_m
  use rint13_m

  implicit none

  integer,  intent(in) :: jatom, LJmax
  real(R8), intent(in) :: bm

  real(R8) :: A(Nrad), B(Nrad), X(Nrad), Y(Nrad), RX
  integer  :: L_index,l1,l2,lj, R_index, i, if1,if2

  do  I=1,JRJ(JATOM)
     RX=R0(JATOM)*exp(DX(JATOM)*(i-1))*BM
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
                 do  I=1,JRJ(JATOM)
                    A(i)=rf1(i,l1,if1)*rj(lj,i)
                    B(i)=rf2(i,l1,if1)*rj(lj,i)
                    X(i)=rf1(i,l2,if2)
                    Y(i)=rf2(i,l2,if2)
                 enddo
                 call RINT13(A,B,X,Y, &
                      ri_mat(r_index,l_index),JATOM)

              end do
           end do
        end do ljloop
     end do
  end do
end subroutine radint
end module radint_m


module     rotate_m; contains
subroutine rotate(VECTOR,ROTMAT,ROTVEC)
  !     ROTATE PERFORMS A ROTATION OF THE VECTOR FROM THE GENERAL
  !     CARTESIAN COORDINATION SYSTEM INTO THE  LOCAL ONE  OF THE
  !     JATOM-TH SPHERE.
  !     THIS SUBROUTINE IS ONLY REQUIRED FOR NONSYMMORPHIC CASES.

  use const, only: R8

  implicit none

  real(R8)    :: VECTOR(3),ROTVEC(3),ROTMAT(3,3)
  intent(in)  :: vector, rotmat
  intent(out) :: rotvec

  integer  :: jcoord, j
  real(R8) :: dotpro

  do JCOORD=1,3
     DOTPRO=0
     do J=1,3
        DOTPRO=DOTPRO + VECTOR(J)*ROTMAT(JCOORD,J)
     end do
     ROTVEC(JCOORD)=DOTPRO
  end do
end subroutine rotate
end module     rotate_m


module     rotdef_m; contains
subroutine rotdef
  !
  !     THE MATRICES  ROTIJ(3,3,INDEX)  TRANSFORM THE POSITION OF AN
  !     ATOM TO ITS CORRESPONDING POSITION OF AN EQUIVALENT ATOM.
  !
  use const,  only: R8
  use struct, only: Nat, mult, iord, iz, pos, trans, transij, rotij, lattic
  use clio,   only: croak
  use util,   only: string

  implicit none

  real(R8)    x(3),x1(3),toler,toler2,one
  INTEGER   i,m,index,index1,jatom,ncount
  DATA TOLER/1.D-7/,ONE/1.D0/
  toler2=1.5d0*toler
  INDEX=0
  NCOUNT=0
  DO JATOM=1,NAT
     INDEX1=INDEX+1
     DO M=1,MULT(JATOM)
        INDEX=INDEX+1
        DO I=1,IORD
           x(1:3)=0.d0
           x=MATMUL(TRANSPOSE(iz(1:3,1:3,i)),pos(1:3,index1))
           x(1:3)=x(1:3)+trans(1:3,i)
           x1(1:3)=MOD(ABS(X(1:3)-POS(1:3,INDEX))+toler,one)-toler
!           WRITE(*,*) 'JATOM,INDEX,I',JATOM,INDEX,I
!           WRITE(*,*) ABS(X1(1:3)),toler
           IF(MAXVAL(ABS(X1)) < TOLER2) THEN
              NCOUNT=NCOUNT+1
              TRANSIJ(1:3,INDEX)=TRANS(1:3,I)
              ROTIJ(1:3,1:3,INDEX)=IZ(1:3,1:3,I)
              GOTO 30
           END IF
           !....check positions for centered lattices
           if(lattic(1:1) == 'B') then
              x1(1:3)=mod(x1(1:3)+0.5d0+toler,one)-toler
              IF(MAXVAL(ABS(X1)) < TOLER2) THEN
                 NCOUNT=NCOUNT+1
                 TRANSIJ(1:3,INDEX)=TRANS(1:3,I)
                 ROTIJ(1:3,1:3,INDEX)=IZ(1:3,1:3,I)
                 GOTO 30
              END IF
           endif
           if(lattic(1:1) == 'F' .or. lattic(1:3) == 'CXY') then
              x1(1:2)=mod(x1(1:2)+0.5d0+toler,one)-toler
              IF(MAXVAL(ABS(X1)) < TOLER2) THEN
                 NCOUNT=NCOUNT+1
                 TRANSIJ(1:3,INDEX)=TRANS(1:3,I)
                 ROTIJ(1:3,1:3,INDEX)=IZ(1:3,1:3,I)
                 GOTO 30
              END IF
              x1(1:2)=mod(x1(1:2)+0.5d0,one)
           endif
           if(lattic(1:1) == 'F' .or. lattic(1:3) == 'CXZ') then
              x1(1)=mod(x1(1)+0.5d0+toler,one)-toler
              x1(3)=mod(x1(3)+0.5d0+toler,one)-toler
              IF(MAXVAL(ABS(X1)) < TOLER2) THEN
                 NCOUNT=NCOUNT+1
                 TRANSIJ(1:3,INDEX)=TRANS(1:3,I)
                 ROTIJ(1:3,1:3,INDEX)=IZ(1:3,1:3,I)
                 GOTO 30
              END IF
              x1(1)=mod(x1(1)+0.5d0,one)
              x1(3)=mod(x1(3)+0.5d0,one)
           endif
           if(lattic(1:1) == 'F' .or. lattic(1:3) == 'CYZ') then
              x1(2:3)=mod(x1(2:3)+0.5d0+toler,one)-toler
              IF(MAXVAL(ABS(X1)) < TOLER2) THEN
                 NCOUNT=NCOUNT+1
                 TRANSIJ(1:3,INDEX)=TRANS(1:3,I)
                 ROTIJ(1:3,1:3,INDEX)=IZ(1:3,1:3,I)
                 GOTO 30
              END IF
           end if
        ENDDO

        call croak("error in ROTDEF: no symmetry operation found &
             &for atoms "//trim(string(jatom))//" and "//trim(string(index))&
             &//", positions "//trim(string(POS(:, index1)))//" and " &
             &//trim(string(POS(:, INDEX))))
30      CONTINUE
     ENDDO
  ENDDO

  if (NCOUNT /= INDEX) call croak('error in ROTDEF: &
       &ROTIJ not defined for all atoms of basis.  NCOUNT = ' &
       &// trim(string(Ncount)))
end subroutine rotdef
end module     rotdef_m


module     latgen_m; contains
subroutine latgen
!
!     LATGEN GENERATES TWO BRAVAIS MATRICES, DEFINES THE VOLUME OF
!     THE UNIT CELL AND CALLS ROTDEF
!     BR1(3,3)  : TRANSFORMS INTEGER RECIPROCAL LATTICE VECTORS AS
!                 GIVEN IN THE VECTORLIST OF LAPW1  ( GENERATED IN
!                 COORS, TRANSFORMED IN BASISO, AND WRITTEN OUT IN
!                 WFTAPE) INTO CARTESIAN SYSTEM
!     BR2(3,3) :  TRANSFORMS A RECIPROCAL LATTICE VECTOR OF A SPE-
!                 CIAL COORDINATE SYSTEM ( IN UNITS OF 2 PI / A )
!                 TO CARTESIAN SYSTEM
  use param,  only: unit_out
  use struct, only: alpha, aa, bb, cc, lattic, pia, vol
  use gener,  only: br1, br2
  use const,  only: R8, TAU
  use clio,   only: croak

  !! procedure includes
  use rotdef_m

  implicit none

  logical  :: ORTHO
  integer  :: i, j
  real(R8) :: rvfac, sinab,sinbc, cosab,cosbc,cosac, wurzel

  real(R8), parameter :: sqrt3=sqrt(3._R8)

  !---------------------------------------------------------------------
  !
  ALPHA(1)=ALPHA(1)*TAU/360
  ALPHA(2)=ALPHA(2)*TAU/360
  ALPHA(3)=ALPHA(3)*TAU/360
  PIA(1)=TAU/AA
  PIA(2)=TAU/BB
  PIA(3)=TAU/CC
  IF(LATTIC(1:1).EQ.'H') GOTO 10
  IF(LATTIC(1:1).EQ.'S') GOTO 20
  IF(LATTIC(1:1).EQ.'P') GOTO 20
  IF(LATTIC(1:1).EQ.'F') GOTO 30
  IF(LATTIC(1:1).EQ.'B') GOTO 40
  IF(LATTIC(1:1).EQ.'C') GOTO 50
  IF(LATTIC(1:1).EQ.'R') GOTO 60

  call croak('error in LATGEN: bad lattice type `'//trim(lattic)//"'")

  !.....HEXAGONAL LATTICE
10 CONTINUE
  BR1(1,1)=2.D0/SQRT3*PIA(1)
  BR1(1,2)=1.D0/SQRT3*PIA(1)
  BR1(1,3)=0.0D0
  BR1(2,1)=0.0D0
  BR1(2,2)=PIA(2)
  BR1(2,3)=0.0D0
  BR1(3,1)=0.0D0
  BR1(3,2)=0.0D0
  BR1(3,3)=PIA(3)
  !
  BR2(1,1)=1.D0
  BR2(1,2)=0.D0
  BR2(1,3)=0.D0
  BR2(2,1)=0.0D0
  BR2(2,2)=1.D0
  BR2(2,3)=0.0D0
  BR2(3,1)=0.0D0
  BR2(3,2)=0.0D0
  BR2(3,3)=1.D0
  !
  RVFAC=2.D0/SQRT(3.D0)
  ORTHO=.FALSE.
  GOTO 100
  !
  !.....RHOMBOHEDRAL CASE
60 BR1(1,1)=1.D0/SQRT(3.D0)*PIA(1)
  BR1(1,2)=1.D0/SQRT(3.D0)*PIA(1)
  BR1(1,3)=-2.d0/sqrt(3.d0)*PIA(1)
  BR1(2,1)=-1.0d0*PIA(2)
  BR1(2,2)=1.0d0*PIA(2)
  BR1(2,3)=0.0d0*PIA(2)
  BR1(3,1)=1.0d0*PIA(3)
  BR1(3,2)=1.0d0*PIA(3)
  BR1(3,3)=1.0d0*PIA(3)
  !
  BR2(1,1)=1.D0
  BR2(1,2)=0.D0
  BR2(1,3)=0.D0
  BR2(2,1)=0.0D0
  BR2(2,2)=1.D0
  BR2(2,3)=0.0D0
  BR2(3,1)=0.0D0
  BR2(3,2)=0.0D0
  BR2(3,3)=1.D0
  !
  RVFAC=6.D0/SQRT(3.D0)
  ORTHO=.FALSE.
  GOTO 100
  !
  !.....PRIMITIVE LATTICE
  !
20 CONTINUE
  SINBC=SIN(ALPHA(1))
  COSAB=COS(ALPHA(3))
  COSAC=COS(ALPHA(2))
  COSBC=COS(ALPHA(1))
  WURZEL=SQRT(SINBC**2-COSAC**2-COSAB**2+2*COSBC*COSAC*COSAB)
  !
  BR1(1,1)= SINBC/WURZEL*PIA(1)
  BR1(1,2)= (-COSAB+COSBC*COSAC)/(SINBC*WURZEL)*PIA(2)
  BR1(1,3)= (COSBC*COSAB-COSAC)/(SINBC*WURZEL)*PIA(3)
  BR1(2,1)= 0.0
  BR1(2,2)= PIA(2)/SINBC
  BR1(2,3)= -PIA(3)*COSBC/SINBC
  BR1(3,1)= 0.0
  BR1(3,2)= 0.0
  BR1(3,3)= PIA(3)
  !
  BR2(1,1)=1.D0
  BR2(1,2)=0.D0
  BR2(1,3)=0.D0
  BR2(2,1)=0.0D0
  BR2(2,2)=1.D0
  BR2(2,3)=0.0D0
  BR2(3,1)=0.0D0
  BR2(3,2)=0.0D0
  BR2(3,3)=1.D0
  !
  RVFAC= 1.d0/WURZEL
  ORTHO=.TRUE.
  if(abs(alpha(1)-TAU/4).gt.0.0001) ortho=.false.
  if(abs(alpha(2)-TAU/4).gt.0.0001) ortho=.false.
  if(abs(alpha(3)-TAU/4).gt.0.0001) ortho=.false.
  !
  GOTO 100
  !
  !.....FC LATTICE
30 CONTINUE
  BR1(1,1)=PIA(1)
  BR1(1,2)=0.0D0
  BR1(1,3)=0.0D0
  BR1(2,1)=0.0D0
  BR1(2,2)=PIA(2)
  BR1(2,3)=0.0D0
  BR1(3,2)=0.0D0
  BR1(3,1)=0.0D0
  BR1(3,3)=PIA(3)
  !
  !     definitions according to column, rows convention for BR2
  !
  BR2(1,1)=-1.D0
  BR2(1,2)= 1.D0
  BR2(1,3)= 1.D0
  BR2(2,1)= 1.D0
  BR2(2,2)=-1.D0
  BR2(2,3)= 1.D0
  BR2(3,1)= 1.D0
  BR2(3,2)= 1.D0
  BR2(3,3)=-1.D0
  !
  RVFAC=4.D0
  ORTHO=.TRUE.
  GOTO 100
  !
  !.....BC LATTICE
40 CONTINUE
  BR1(1,1)=PIA(1)
  BR1(1,2)=0.0D0
  BR1(1,3)=0.0D0
  BR1(2,1)=0.0D0
  BR1(2,2)=PIA(2)
  BR1(2,3)=0.0D0
  BR1(3,1)=0.0D0
  BR1(3,2)=0.0D0
  BR1(3,3)=PIA(3)
  !
  BR2(1,1)= 0.0D0
  BR2(1,2)= 1.D0
  BR2(1,3)= 1.D0
  BR2(2,1)= 1.D0
  BR2(2,2)= 0.0D0
  BR2(2,3)= 1.D0
  BR2(3,1)= 1.D0
  BR2(3,2)= 1.D0
  BR2(3,3)= 0.0D0
  !
  RVFAC=2.D0
  ORTHO=.TRUE.
  GOTO 100
  !
50 CONTINUE
  IF(LATTIC(2:3).EQ.'XZ') GOTO 51
  IF(LATTIC(2:3).EQ.'YZ') GOTO 52
  !.....CXY LATTICE
  BR1(1,1)=PIA(1)
  BR1(1,2)=0.0D0
  BR1(1,3)=0.0D0
  BR1(2,1)=0.0D0
  BR1(2,2)=PIA(2)
  BR1(2,3)=0.0D0
  BR1(3,1)=0.0D0
  BR1(3,2)=0.0D0
  BR1(3,3)=PIA(3)
  !
  BR2(1,1)= 1.D0
  BR2(1,2)= 1.D0
  BR2(1,3)= 0.0D0
  BR2(2,1)=-1.D0
  BR2(2,2)= 1.D0
  BR2(2,3)= 0.0D0
  BR2(3,1)= 0.0D0
  BR2(3,2)= 0.0D0
  BR2(3,3)= 1.D0
  !
  RVFAC=2.D0
  ORTHO=.TRUE.
  GOTO 100
  !
  !.....CXZ CASE (CXZ LATTICE BUILD UP)
51 CONTINUE
  !.....CXZ ORTHOROMBIC CASE
  if(abs(ALPHA(3)-TAU/4).lt.0.0001) then
     BR1(1,1)=PIA(1)
     BR1(1,2)=0.0D0
     BR1(1,3)=0.0D0
     BR1(2,1)=0.0D0
     BR1(2,2)=PIA(2)
     BR1(2,3)=0.0D0
     BR1(3,1)=0.0D0
     BR1(3,2)=0.0D0
     BR1(3,3)=PIA(3)
     !
     BR2(1,1)= 1.D0
     BR2(1,2)= 0.0
     BR2(1,3)= 1.D0
     BR2(2,1)= 0.0
     BR2(2,2)= 1.D0
     BR2(2,3)= 0.0
     BR2(3,1)=-1.D0
     BR2(3,2)= 0.0
     BR2(3,3)= 1.D0
     !
     RVFAC=2.0
     ORTHO=.TRUE.
     GOTO 100
  ELSE
     !.....CXZ MONOCLINIC CASE
     write(*,*) '  gamma not equal 90'
     SINAB=SIN(ALPHA(3))
     COSAB=COS(ALPHA(3))
     !
     BR1(1,1)= PIA(1)/SINAB
     BR1(1,2)= -PIA(2)*COSAB/SINAB
     BR1(1,3)= 0.0
     BR1(2,1)= 0.0
     BR1(2,2)= PIA(2)
     BR1(2,3)= 0.0
     BR1(3,1)= 0.0
     BR1(3,2)= 0.0
     BR1(3,3)= PIA(3)
     !
     BR2(1,1)= 1.D0
     BR2(1,2)= 0.D0
     BR2(1,3)= 1.D0
     BR2(2,1)= 0.0
     BR2(2,2)= 1.D0
     BR2(2,3)= 0.0
     BR2(3,1)=-1.D0
     BR2(3,2)= 0.0
     BR2(3,3)= 1.D0
     !
     RVFAC=2.0/SINAB
     ORTHO=.FALSE.
     GOTO 100
  ENDIF
  !
  !.....CYZ CASE (CYZ LATTICE BUILD UP)
52 CONTINUE
  BR1(1,1)=PIA(1)
  BR1(1,2)=0.0D0
  BR1(1,3)=0.0D0
  BR1(2,1)=0.0D0
  BR1(2,2)=PIA(2)
  BR1(2,3)=0.0D0
  BR1(3,1)=0.0D0
  BR1(3,2)=0.0D0
  BR1(3,3)=PIA(3)
  !
  BR2(1,1)= 1.D0
  BR2(1,2)= 0.0
  BR2(1,3)= 0.0
  BR2(2,1)= 0.0
  BR2(2,2)= 1.D0
  BR2(2,3)= 1.D0
  BR2(3,1)= 0.0
  BR2(3,2)=-1.D0
  BR2(3,3)= 1.D0
  !
  RVFAC=2.0
  ORTHO=.TRUE.
  GOTO 100
  !
  !.....DEFINE VOLUME OF UNIT CELL
100 CONTINUE
  write(unit_out,*)' BR1,  BR2'
  do i=1,3
     write(unit_out,654)(br1(i,j),j=1,3),(br2(i,j),j=1,3)
654  format(3f10.5,3x,3f10.5)
  enddo
  VOL=AA*BB*CC/RVFAC
  !
  !.....DEFINE ROTATION MATRICES IN NONSYMMORPHIC CASE
  call ROTDEF
end subroutine latgen
end module     latgen_m


!!/---
!! Local Variables:
!! mode: f90
!! End:
!!\---
!!
!! Time-stamp: <2016-07-25 16:08:21 assman@faepop71.tu-graz.ac.at>
