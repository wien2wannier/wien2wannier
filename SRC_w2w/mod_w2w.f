!!! wien2wannier/SRC_w2w/mod_w2w.f
!!!
!!!    Collection of routines, variables, and constants shared between
!!!    w2w and wplot.

!!/=== Defined here: =============================
!!
!! param:   wien2wannier_version, Lmax2, LOmax, Nloat, Nrad, Nrf,
!!          unit_def, unit_in, unit_out, unit_struct, unit_vector,
!!          unit_vsp
!!
!! PS1:     DEP(5), DEQ(5), DB, DVC, DSAL, DK, DM
!!
!! uhelp:   A(Nrad), B(Nrad)
!!
!! Procedure modules:
!!
!!    dergl_m, diracout_m, dvbes1_m, inouh_m, inth_m, sphbes_m, Ylm_m
!!
!!\===============================================

module param
  use const, only: R8, C16
  implicit none
  private :: R8, C16

  public

  character(*), parameter, private :: &
       rev_str = "$version: v1.0.0-213-g12b1846$"
  character(*), parameter, public  :: &
       wien2wannier_version = rev_str(11 : len (rev_str)-1)

  integer, parameter :: unit_def=1, unit_in=5, unit_out=6, unit_vector=10
  integer, parameter :: unit_vsp=18, unit_struct=20
  
  integer, parameter :: Lmax2=5, LOmax=3, Nloat=3, Nrad=881, Nrf=4
end module param

module PS1
  use const, only: R8
  implicit none
  private; save

  ! DEP,DEQ DERIVEES DE DP ET DQ   DB=ENERGIE/DVC    DVC VITESSE DE LA
  ! LUMIERE EN U.A.   DSAL=2.*DVC   DK NOMBRE QUANTIQUE KAPPA
  ! DM=PAS EXPONENTIEL/720., DKOEF=1./720.
  real(R8), public :: DEP(5), DEQ(5), DB, DVC, DSAL, DK, DM
end module PS1

module uhelp
  use param, only: Nrad
  use const, only: R8

  implicit none
  private; save

  real(R8), public :: A(Nrad), B(Nrad)
end module uhelp


!---------------------------  Procedure modules  ---------------------------
module     dergl_m; contains
!rschmid
!     Calculate the first derivate of f_l.
subroutine dergl(stru, jatom, g_l, der1)
  use param,     only: Nrad
  use const,     only: R8
  use structmod, only: struct_t

  implicit none

  type(struct_t), intent(in)  :: stru
  integer,        intent(in)  :: jatom
  real(R8),       intent(in)  :: g_l(nrad)
  real(R8),       intent(out) :: der1(nrad)

  real(R8) :: f_l(nrad), dfldi(nrad), r(nrad), rad
  integer  :: j,i

  do i = 1,stru%Npt(jatom)
     r(i) = stru%R0(jatom) * exp(stru%dx(jatom) * (i-1))
     f_l(i) = g_l(i)/r(i)
  enddo

  !rschmid
  !  Use unsymmetric 6-point formulae for the first derivative at the
  !  first 3 mesh points.
  !rschmid
  dfldi(1) =  ( - 137*f_l(1) + 300*f_l(2)   &
                - 300*f_l(3) + 200*f_l(4)   &
                -  75*f_l(5) +  12*f_l(6) ) &
                 / 60
  dfldi(2) = ( -  12*f_l(1) -  65*f_l(2)    &
               + 120*f_l(3) -  60*f_l(4)    &
               +  20*f_l(5) -   3*f_l(6) )  &
                 / 60
  dfldi(3) = (     3*f_l(1) -  30*f_l(2)    &
               -  20*f_l(3) +  60*f_l(4)    &
               -  15*f_l(5) +   2*f_l(6) )  &
                 / 60

  RAD=R(1)
  DER1(1)=dfldi(1)/RAD/stru%dx(jatom)
  RAD=R(2)
  DER1(2)=dfldi(2)/RAD/stru%dx(jatom)
  RAD=R(3)
  DER1(3)=dfldi(3)/RAD/stru%dx(jatom)
  !rschmid
  !  Use symmetric 7-point formula to generate the first derivative at
  !  all intermediate mesh points.
  !rschmid

  do J=4,stru%Npt(jatom)-3
     dfldi(J) = (        f_l(J+3) - f_l(J-3)    &
                 -  9 * (f_l(J+2) - f_l(J-2))   &
                 + 45 * (f_l(J+1) - f_l(J-1)) ) &
                 / 60
     RAD=R(J)
     DER1(J)=dfldi(J)/RAD/stru%dx(jatom)
  end do
!rschmid
!  Use unsymmetric 6-point formulae for the second derivative at the
!  last 3 mesh points.
!rschmid
  dfldi(stru%Npt(jatom)-2) = - (    3*f_l(stru%Npt(jatom)  )   &
                                -  30*f_l(stru%Npt(jatom)-1)   &
                                -  20*f_l(stru%Npt(jatom)-2)   &
                                +  60*f_l(stru%Npt(jatom)-3)   &
                                -  15*f_l(stru%Npt(jatom)-4)   &
                                +   2*f_l(stru%Npt(jatom)-5) ) &
                                / 60
  RAD=R(stru%Npt(jatom)-2)
  DER1(stru%Npt(jatom)-2)=dfldi(stru%Npt(jatom)-2)/RAD/stru%dx(jatom)

  dfldi(stru%Npt(jatom)-1) = - (-  12*f_l(stru%Npt(jatom)  )   &
                                -  65*f_l(stru%Npt(jatom)-1)   &
                                + 120*f_l(stru%Npt(jatom)-2)   &
                                -  60*f_l(stru%Npt(jatom)-3)   &
                                +  20*f_l(stru%Npt(jatom)-4)   &
                                -   3*f_l(stru%Npt(jatom)-5) ) &
                                / 60
  RAD=R(stru%Npt(jatom)-1)
  DER1(stru%Npt(jatom)-1)=dfldi(stru%Npt(jatom)-1)/RAD/stru%dx(jatom)

  dfldi(stru%Npt(jatom)  ) = - (- 137*f_l(stru%Npt(jatom)  )   &
                                + 300*f_l(stru%Npt(jatom)-1)   &
                                - 300*f_l(stru%Npt(jatom)-2)   &
                                + 200*f_l(stru%Npt(jatom)-3)   &
                                -  75*f_l(stru%Npt(jatom)-4)   &
                                +  12*f_l(stru%Npt(jatom)-5) ) &
                                / 60

  RAD=R(stru%Npt(jatom))
  DER1(stru%Npt(jatom)) = dfldi(stru%Npt(jatom))/RAD/stru%dx(jatom)
end subroutine dergl
end module     dergl_m


module     inouh_m; contains
subroutine inouh(dp, dq, dr, dq1, dfl, dv, z, test, nuc)
!!! valeurs initiales pour l integration vers l exrerieur
!!!    dp    grande composante
!!!    dq    petite composante
!!!    dr    bloc des points
!!!    dq1   pente a l origine de dp ou dq
!!!    dfl   puissance du premier terme
!!!    dv    potentiel au premier point
!!!    z     numero atomique
!!!    test  test de precision
!!!    nuc   noyaux de dimensions finies si nuc non nul

  use param, only: Nrad
  use const, only: R8
  use PS1,   only: dep, deq, db, dvc, dsal, dk ! we define out own dm

  implicit none

  integer  :: nuc
  real(R8) :: dp(nrad), dq(nrad), dr(nrad), dq1, dfl, dv, z, test

  intent(in)  :: nuc, dr, dq1, dfl, dv, z, test
  intent(out) :: dp, dq

  integer  :: i,j, m
  real(R8) :: dval, dm, deva1,deva2,deva3, dbe, dsum, dpr,dqr

!!! dep,deq   derivees de dp et dq
!!! db    =   energie/dvc
!!! dvc       vitesse de la lumiere en u.a.
!!! dsal  =   2.*dvc
!!! dk        nombre quantique kappa
!!! dm        pas exponentiel/720.

      do i=1,10
        dp(i)=0.
        dq(i)=0.
      enddo

      if (nuc.le.0) then
        dval=z/dvc
        deva1=-dval
        deva2=dv/dvc+dval/dr(1)-db
        deva3=0.
        if (dk.le.0) then
          dbe=(dk-dfl)/dval
        else
          dbe=dval/(dk+dfl)
        endif
        dq(10)=dq1
        dp(10)=dbe*dq1

      else
        dval=dv+z*(3-dr(1)*dr(1)/(dr(nuc)*dr(nuc)))/(dr(nuc)+dr(nuc))
        deva1=0
        deva2=(dval-3*z/(dr(nuc)+dr(nuc)))/dvc-db
        deva3=z/(dr(nuc)*dr(nuc)*dr(nuc)*dsal)
        if (dk.le.0.d0) then
          dp(10)=dq1
        else
          dq(10)=dq1
        endif
      endif
      do i=1,5
        dp(i)=dp(10)
        dq(i)=dq(10)
        dep(i)=dp(i)*dfl
        deq(i)=dq(i)*dfl
      enddo
      m=1
 41   dm=m+dfl
      dsum=dm*dm-dk*dk+deva1*deva1
      dqr=(dsal-deva2)*dq(m+9)-deva3*dq(m+7)
      dpr=deva2*dp(m+9)+deva3*dp(m+7)
      dval=((dm-dk)*dqr-deva1*dpr)/dsum
      dsum=((dm+dk)*dpr+deva1*dqr)/dsum
      j=-1
      do i=1,5
        dpr=dr(i)**m
        dqr=dsum*dpr
        dpr=dval*dpr
        if (m.ne.1) then
          if (abs(dpr/dp(i)).le.test.and.abs(dqr/dq(i)).le.test) j=1
        endif
        dp(i)=dp(i)+dpr
        dq(i)=dq(i)+dqr
        dep(i)=dep(i)+dpr*dm
        deq(i)=deq(i)+dqr*dm
      enddo
      if (j.ne.1) then
        dp(m+10)=dval
        dq(m+10)=dsum
        m=m+1
        if (m <= 20) goto 41
      endif
      return
    end subroutine inouh
end module inouh_m


module     inth_m; contains
! INTEGRATION PAR LA METHODE DE ADAMS A 5 POINTS DE LA GRANDE COMPOSANTE
! DP ET DE LA PETITE COMPOSANTE DQ AU POINT DR,DV ETANT LE POTENTIEL EN
! CE POINT
subroutine inth(DP,DQ,DV,DR)
  use const, only: R8
  use PS1,   only: DEP, DEQ, DB, DVC, DSAL, DK, DM

  implicit none

  real(R8), intent(inout) :: dp, dq
  real(R8), intent(in)    :: dv, dr

!!! DEP,DEQ DERIVEES DE DP ET DQ
!!! DB    = ENERGIE/DVC
!!! DVC     VITESSE DE LA LUMIERE EN U.A.
!!! DSAL  = 2.*DVC
!!! DK      NOMBRE QUANTIQUE KAPPA
!!! DM    = PAS EXPONENTIEL/720.
!!! DKOEF1= 405./502. (!)
!!! DKOEF2= 27. /502.

  ! dkoef1 is inferred from old numerical literal
  real(R8), parameter :: dkoef1 = 475 / 502._R8
  real(R8), parameter :: dkoef2 =  27 / 502._R8

  real(R8) :: DPR, DQR, dsum
  integer  :: i

  DPR=DP+DM*((251.*DEP(1)+2616.*DEP(3)+1901.*DEP(5))-(1274.*DEP(2)+2774.*DEP(4)))
  DQR=DQ+DM*((251.*DEQ(1)+2616.*DEQ(3)+1901.*DEQ(5))-(1274.*DEQ(2)+2774.*DEQ(4)))

  do I=2,5
     DEP(I-1)=DEP(I)
     DEQ(I-1)=DEQ(I)
  end do

  DSUM=(DB-DV/DVC)*DR

  DEP(5)=-DK*DPR + (DSAL*DR+DSUM)*DQR
  DEQ(5)= DK*DQR -  DSUM*DPR

  DP=DP+DM*((106.*DEP(2)+646.*DEP(4)+251.*DEP(5))-(19.*DEP(1)+264.*DEP(3)))
  DQ=DQ+DM*((106.*DEQ(2)+646.*DEQ(4)+251.*DEQ(5))-(19.*DEQ(1)+264.*DEQ(3)))
  DP = DKOEF1 * DP + DKOEF2 * DPR
  DQ = DKOEF1 * DQ + DKOEF2 * DQR
  DEP(5) = -DK * DP + (DSAL*DR+DSUM)*DQ
  DEQ(5) =  DK * DQ - DSUM * DP
end subroutine inth
end module     inth_m


module     diracout_m; contains
!!! Integration of Dirac equation
subroutine diracout(stru, jatom, V, eh, nqk, val, slo, nodes)
  use param,     only: unit_out, Nrad
  use const,     only: clight, R8
  ! dp    = large component of the solution of the dirac equation
  ! dq    = small component of the solution
  use uhelp,     only: dp => A, dq => B
  use PS1,       only: dep, deq, db, dvc, dsal, dk, dm
  use structmod, only: struct_t

  !! procedure includes
  use inth_m
  use inouh_m

  implicit none

  !  Input:
  !    stru   structure data type
  !    V      rad.sym. potential in Hartree, V = potential*r
  !    Nmax   number of radial meshpoints
  !    EH     energy in hartree
  !    Nqk    relativistic quantum number kappa
  !    Z      charge of nucleus
  type(struct_t), intent(in) :: stru
  real(R8),       intent(in) :: V(Nrad), EH
  integer,        intent(in) :: jatom, Nqk

  !  Output:
  !    val,slo:  Wellenfunktion und Steigung am Kugelrand
  !    nodes:    nomber of nodes
  real(R8), intent(out) :: val
  integer,  intent(out) :: nodes

  real(R8), parameter :: dkoef=1 / 720._R8, test=1.e-8_R8
  real(R8), parameter :: atom_mass(103) = (/ &
       &   1.0_R8,   4.0_R8,   6.9_R8,   9.0_R8,  10.8_R8,  12.0_R8, &
       &  14.0_R8,  16.0_R8,  19.0_R8,  20.2_R8,  23.0_R8,  24.3_R8, &
       &  27.0_R8,  28.1_R8,  31.0_R8,  32.0_R8,  35.4_R8,  40.0_R8, &
       &  39.1_R8,  40.0_R8,  45.0_R8,  47.9_R8,  50.9_R8,  52.0_R8, &
       &  54.9_R8,  55.8_R8,  58.9_R8,  58.7_R8,  63.5_R8,  65.4_R8, &
       &  69.7_R8,  72.6_R8,  74.9_R8,  79.0_R8,  79.9_R8,  83.8_R8, &
       &  85.5_R8,  87.6_R8,  88.9_R8,  91.2_R8,  92.9_R8,  95.9_R8, &
       &  98.0_R8, 101.1_R8, 102.9_R8, 106.4_R8, 107.9_R8, 112.4_R8, &
       & 114.8_R8, 118.7_R8, 121.8_R8, 127.6_R8, 126.9_R8, 131.3_R8, &
       & 132.9_R8, 137.3_R8, 138.9_R8, 140.1_R8, 140.9_R8, 144.2_R8, &
       & 145.0_R8, 150.4_R8, 152.0_R8, 157.3_R8, 158.9_R8, 162.5_R8, &
       & 164.9_R8, 167.3_R8, 168.9_R8, 173.0_R8, 175.0_R8, 178.5_R8, &
       & 180.9_R8, 183.8_R8, 186.2_R8, 190.2_R8, 192.2_R8, 195.1_R8, &
       & 197.0_R8, 200.6_R8, 204.4_R8, 207.2_R8, 209.0_R8, 209.0_R8, &
       & 210.0_R8, 222.0_R8, 223.0_R8, 226.0_R8, 227.0_R8, 232.0_R8, &
       & 231.0_R8, 238.0_R8, 237.0_R8, 244.0_R8, 243.0_R8, 247.0_R8, &
       & 247.0_R8, 251.0_R8, 252.0_R8, 257.0_R8, 258.0_R8, 259.0_R8, &
       & 262.0_R8 /)

  !      DR   =    radial mesh
  real(R8) :: dv(nrad),  dr(NRAD)
  real(R8) :: d1, dfl, dq1, Rnuc, dval, slo
  integer  :: i, nuc

  !rschmid
  !   Set up radial mesh.
  !rschmid
  do i=1,stru%Npt(jatom)
     DR(i)=stru%R0(jatom) * (exp(stru%dx(jatom) * (i-1.d0)))
  enddo

  if (stru%rel) then
     dvc = clight
  else
     dvc = 1e10_R8
  endif
  dsal = 2*dvc
  db = eh/dvc
  dk = nqk
  dm = stru%dx(jatom) * dkoef

  do i=1,stru%Npt(jatom)
     dv(i) = v(i)/dr(i)
  enddo
  !rschmid
  !  Behavior of the solution at the origin
  !rschmid

  !jk   finite size of the nucleus

  ! Shouldn't we use nint() here instead of int()?  OTOH, all of
  ! Wien2k does it this way …
  Rnuc = 2.2677e-05_R8 * atom_mass(int(stru%Z(jatom)))**(1/3._R8)
  write(unit_out,*) 'amass, r0:', atom_mass(int(stru%Z(jatom))), Rnuc
  do i = 1, stru%Npt(jatom)
     d1 = stru%R0(jatom) * exp(stru%dx(jatom) * (i-1.d0))
     if (d1 >= Rnuc) exit
  end do
  nuc=I
  write(unit_out,*)'nuc=',nuc
  if (nuc <= 0) then
     dfl = sqrt(nqk*nqk - stru%Z(jatom)**2 / dvc**2)
  else
     dfl=nqk*nqk
     do i=1,nuc
        dv(i) = dv(i) + stru%Z(jatom)/dr(i) + &
             &  stru%Z(jatom) * ((dr(i)/dr(nuc))**2 - 3) / (2*dr(nuc))
     end do
  end if
  dq1 = nqk/iabs(nqk)

  !rschmid
  !  Determine expansion of the potential at the origin.
  !rschmid
  call inouh(dp, dq, dr, dq1, dfl, dv(1), stru%Z(jatom), TEST, nuc)

  !rschmid
  !  Set up boundary conditions using a small r expansion.
  !rschmid

  nodes = 0
  do i=1,5
     dval=dr(i)**dfl
     if (i /= 1) then
        if (dp(i-1) /= 0) then
           if ((dp(i)/dp(i-1)) <= 0) then
              nodes=nodes+1
           endif
        endif
     endif
     dp(i) = dp(i)*dval
     dq(i) = dq(i)*dval
     dep(i)=dep(i)*dval
     deq(i)=deq(i)*dval
  enddo

  !rschmid
  !    Perform outward integration of dirac equation
  !rschmid

  do i = unit_out, stru%Npt(jatom)
     dp(i) = dp(i-1)
     dq(i) = dq(i-1)
     call inth (dp(i),dq(i),dv(i),dr(i))
     if (dp(i-1) /= 0) then
        if ((dp(i)/dp(i-1)) > 0) then
           nodes=nodes+1
        endif
     endif
  enddo


  val = dp(stru%Npt(jatom)) / dr(stru%Npt(jatom))
  slo = dep(5) / (stru%dx(jatom) * dr(stru%Npt(jatom))) / dvc*2
  slo = (slo-val) / dr(stru%Npt(jatom))
end subroutine diracout
end module     diracout_m


module     dvbes1_m; contains
!-----X----X----X----X----X----X----X----X----X----X----X----X----X----X
!-----X CALCULATE THE DERIVATIVES OF THE BESSEL FUNCTIONS.   X----X----X
!-----X   DJ=DFJ/DX WHERE X=SM*RI                                 X----X
!-----X                    D.D.KOELLING                      X----X----X
!-----X----X----X----X----X----X----X----X----X----X----X----X----X----X
subroutine dvbes1(FJ,DJ,SM,NT)
  use const, only: DPk

  implicit none

  real(DPk), intent(in)  :: FJ(:), SM
  real(DPk), intent(out) :: DJ(:)
  integer,   intent(in)  :: NT

  real(DPk), parameter :: ZUP = 1e-5_DPk

  integer   :: l, lm
  real(DPk) :: Q2, Q3

  if(SM <= ZUP) then
     DJ(1) = 0
     DJ(2) = 1/3._DPk
     do L=3,NT
        DJ(L) = 0
     end do
  else
     Q2=-1/SM
     Q3=Q2
     DJ(1)=-FJ(2)
     LM=1
     do L=2,NT
        Q3=Q3+Q2
        DJ(L)=FJ(LM)+Q3*FJ(L)
        LM=LM+1
     end do
  end if
end subroutine dvbes1
end module     dvbes1_m


module     sphbes_m; contains
!***  VERSION III-UPPER LIMIT OBTAINED FROM THE EXPRESSIONS OF
!***             CORBATO AND URETSKY USING A LIMIT OF E SUB M OF
!***             2**-30 WHICH IS APPROXIMATELY 1.E-9
!***            SUMMATION PROPERTY USED TO NORMALIZE RESULTS.
!***  ADDITIONAL FACTOR ADDED TO STARTING VALUE
!***  N IS THE MAXIMUM L TO BE CALCULATED
!***  X IS THE ARGUMENT
!***  FJ IS THE ARRAY THAT THE SPHERICAL BESSEL FUNCTIONS ARE TO BE
!***  PLACED IN.
!*****  MODIFIED TO NOT REQUIRE THE WORKING SPACE.
!*****        29 MAY,1968
subroutine sphbes(N,X,FJ)
  use const, only: R8
  use clio,  only: croak
  use util,  only: string

  implicit none

  integer,  intent(in)  :: N
  real(R8), intent(in)  :: X
  real(R8), intent(out) :: FJ(:)

  real(R8), parameter :: XLIM = 0.1_R8
  real(R8), parameter :: HF   = 0.5_R8
  real(R8), parameter :: TNHF = 10.5_R8
  real(R8), parameter :: T25  = 1.0e25_R8
  real(R8), parameter :: TN25 = 1.0e-25_R8
  real(R8), parameter :: TN50 = 1.0e-50_R8

  real(R8) :: hfxsq, xl, twm, ta, xlp, cufac, ffo, ffn, xi, fm, sdr, ffp, ser
  integer  :: ns, m, mm, j, jj

  if (N < 0) call croak('error in SPHBES: N should not be negative - ' &
       // trim(string(N)))
  if (X < 0) call croak('error in SPHBES: x should not be negative - ' &
       // trim(string(X)))
  if (X > XLIM) GO TO 25
  HFXSQ=HF*X*X
  XL=1
  TWM=1
  M=0
11 M=M+1
  TA=XL
  TWM=TWM+2
  XL=XL/TWM
  TA=TA-XL*HFXSQ
  XLP=XL/(TWM+2)
  FJ(M)=TA+HF*XLP*HFXSQ*HFXSQ
  XL=XL*X
  IF (M.LE.N)  GO TO 11
  RETURN
25 CUFAC=4.2_R8
  IF (X.LT.(N-2)) CUFAC=TNHF/(N+HF-X)
  NS=N+5+int(X*CUFAC)
  !*******************  ADD ADDITIONAL FACTOR  ***************************
  NS=NS + int(15/(1 + sqrt(X)))
  !***********************************************************************
  CONTINUE
  FFO=0
  FFN=TN25
  M=NS-1
  XI=1/X
  FM=(M+M)+1
  SDR=FM*TN50
314 FFP=FM*XI*FFN-FFO
  IF (ABS(FFP).LT.T25) GO TO 315
  SDR=SDR*TN50
  FFP=FFP*TN25
  FFN=FFN*TN25
315 SDR=SDR + (FM-2)*FFP*FFP
  FFO=FFN
  FFN=FFP
  IF (M.LE.N) GO TO 316
  M=M-1
  FM=FM-2
  GO TO 314
316 FJ(M)=FFN
  FJ(M+1)=FFO
  GO TO 33
32 FJ(M)=FM*XI*FJ(M+1)-FJ(M+2)
  IF(ABS(FJ(M)).GE.T25) GO TO 56
  SDR=SDR + (FM-2)*FJ(M)*FJ(M)
  IF (M.LE.1) GO TO 34
33 M = M-1
  FM=FM-2
  GO TO 32
34 SER=1/SQRT(SDR)
  MM = N+1
  do M=1,MM
     FJ(M)=FJ(M)*SER
  end do
  return
56 JJ= M+1
  NS=N+1
  do J = JJ,NS
     FJ(J)=FJ(J)*TN25
  end do
  SDR=SDR*TN50
  GO TO 32
end subroutine sphbes
end module     sphbes_m


module     Ylm_m; contains
!!! Calculate spherical harmonics
!!!
!!!   The spherical harmonics (Condon and Shortley convention)
!!!     Y(0,0),Y(1,-1),Y(1,0),Y(1,1),Y(2,-2) ... Y(LMAX,LMAX)
!!!   for vector V (given in Cartesian coordinates)
!!!   are calculated. In the Condon Shortley convention the
!!!   spherical harmonics are defined as
!!!                     +------+
!!!                m    |   1     m              im(Phi)
!!!   Y(l,m) = (-1)  -+ | -----  P (cos(Theta)) e
!!!                    \| 2(Pi)   l
!!!          m
!!!   where P (cos(Theta)) is the normalized Associated Legendre
!!!          l
!!!   function. Thus,
!!!                                  m      *
!!!                    Y(l,-m) = (-1) Y(l,m)
!!!
!!! Usage:
!!!   DOUBLE PRECISION V(3), Y(5*5)
!!!   V(1) = ...
!!!   V(2) = ...
!!!   V(3) = ...
!!!   CALL YLM(V,4,Y)
!!!
!!! Arguments:
!!!     V      - DOUBLE PRECISION vector, dimension 3        (input)
!!!              Must be given in Cartesian coordinates.
!!!              Conversion of V to polar coordinates gives the
!!!              angles Theta and Phi necessary for the calculation
!!!              of the spherical harmonics.
!!!     LMAX   - INTEGER value                               (input)
!!!              upper bound of L for which spherical harmonics
!!!              will be calculated
!!!              constraint:
!!!                 LMAX .GE. 0 (not checked)
!!!     Y      - COMPLEX*16 array, dimension (LMAX+1)**2    (output)
!!!              contains the calculated spherical harmonics
!!!              Y(1)                   for L .EQ. 0 (M = 0)
!!!              Y(2), ..., Y(4)        for L .EQ. 1 (M = -1, 0, 1)
!!!              ...
!!!              Y(LMAX*LMAX+1), ..., Y((LMAX+1)*(LMAX+1))
!!!                                     for L .EQ. LMAX
!!!                                        (M = -L,...,L)
!!!              constraint:
!!!                 Dimension of Y .GE. (LMAX+1)**2 (not checked)
!!!
!!! Method:
!!!
!!!   The basic algorithm used to calculate the spherical
!!!   harmonics for vector V is as follows:
!!!
!!!   Y(0,0)
!!!   Y(1,0)
!!!   Y(1,1)
!!!   Y(1,-1) = -Y(1,1)
!!!   DO L = 2, LMAX
!!!      Y(L,L)   = f(Y(L-1,L-1)) ... Formula 1
!!!      Y(L,L-1) = f(Y(L-1,L-1)) ... Formula 2
!!!      DO M = L-2, 0, -1
!!!         Y(L,M) = f(Y(L-1,M),Y(L-2,M)) ... Formula 2
!!!         Y(L,-M)= (-1)**M*Y(L,M)
!!!      ENDDO
!!!   ENDDO
!!!
!!!   In the following the necessary recursion formulas and
!!!   starting values are given:
!!!
!!!  Start:
!!!                  +------+
!!!                  |   1
!!!     Y(0,0) =  -+ | -----
!!!                 \| 4(Pi)
!!!
!!!                             +------+
!!!                             |   3
!!!     Y(1,0) =  cos(Theta) -+ | -----
!!!                            \| 4(Pi)
!!!
!!!                               +------+
!!!                               |   3    i(Phi)
!!!     Y(1,1) =  - sin(Theta) -+ | ----- e
!!!                              \| 8(Pi)
!!!
!!!  Formula 1:
!!!
!!!     Y(l,l) =
!!!                     +--------+
!!!                     | (2l+1)   i(Phi)
!!!      -sin(Theta) -+ | ------  e       Y(l-1,l-1)
!!!                    \|   2l
!!!
!!!  Formula 2:
!!!                            +---------------+
!!!                            |  (2l-1)(2l+1)
!!!     Y(l,m) = cos(Theta) -+ | -------------- Y(l-1,m)  -
!!!                           \|   (l-m)(l+m)
!!!
!!!                              +--------------------+
!!!                              |(l-1+m)(l-1-m)(2l+1)
!!!                        -  -+ |-------------------- Y(l-2,m)
!!!                             \|  (2l-3)(l-m)(l+m)
!!!
!!!  Formula 3: (not used in the algorithm because of the division
!!!              by sin(Theta) which may be zero)
!!!
!!!                              +--------------+
!!!                cos(Theta)    |  4(m+1)(m+1)   -i(Phi)
!!!     Y(l,m) = - ---------- -+ | ------------  e       Y(l,m+1) -
!!!                sin(Theta)   \| (l+m+1)(l-m)
!!!
!!!                              +--------------+
!!!                              |(l-m-1)(l+m+2)  -2i(Phi)
!!!                        -  -+ |-------------- e        Y(l,m+2)
!!!                             \| (l-m)(l+m+1)
subroutine Ylm(V, lmax, Y)
  use const, only: R8, C16, TAU

  implicit none

  integer,      intent(in)  :: LMAX
  real(R8),     intent(in)  :: V(3)
  complex(C16), intent(out) :: Y(:)

  integer  :: I2L, I4L2, INDEX, INDEX2, L, M, MSIGN
  real(R8) :: A, B, C, AB, ABC, ABMAX, ABCMAX, D4LL1C, D2L13
  real(R8) :: COSTH, SINTH, COSPH, SINPH, TEMP1, TEMP2, TEMP3
  real(R8) :: YLLR, YLL1R, YL1L1R, YLMR, YLLI, YLL1I, YL1L1I, YLMI

! Y(0,0)
  YLLR = 1/sqrt(2*TAU)
  YLLI = 0
  Y(1) = cmplx(YLLR, YLLI, R8)

! continue only if spherical harmonics for (L .GT. 0) are desired
  if (LMAX <= 0) return

! calculate sin(Phi), cos(Phi), sin(Theta), cos(Theta)
! Theta, Phi ... polar angles of vector V
  ABMAX  = max(abs(V(1)),abs(V(2)))
  if (ABMAX > 0) then
     A = V(1)/ABMAX
     B = V(2)/ABMAX
     AB = sqrt(A*A+B*B)
     COSPH = A/AB
     SINPH = B/AB
  else
     COSPH = 1
     SINPH = 0
  endif
  ABCMAX = max(ABMAX,abs(V(3)))
  if (ABCMAX > 0) then
     A = V(1)/ABCMAX
     B = V(2)/ABCMAX
     C = V(3)/ABCMAX
     AB = A*A + B*B
     ABC = sqrt(AB + C*C)
     COSTH = C/ABC
     SINTH = sqrt(AB)/ABC
  else
     COSTH = 1
     SINTH = 0
  endif

! Y(1,0)
  Y(3) = cmplx(sqrt(3._R8)*YLLR*COSTH, 0, R8)

! Y(1,1) ( = -DCONJG(Y(1,-1)))
  TEMP1 = -sqrt(1.5D+0)*YLLR*SINTH
  Y(4) = DCMPLX(TEMP1*COSPH,TEMP1*SINPH)
  Y(2) = -DCONJG(Y(4))

  do L = 2, LMAX
     INDEX  = L*L+1
     INDEX2 = INDEX + 2*L
     MSIGN  = 1 - 2*mod(L,2)

!    YLL = Y(L,L) = f(Y(L-1,L-1)) ... Formula 1
     YL1L1R = dble(Y(INDEX-1))
     YL1L1I = aimag(Y(INDEX-1))
     TEMP1 = -sqrt(dble(2*L+1)/dble(2*L))*SINTH
     YLLR = TEMP1*(COSPH*YL1L1R - SINPH*YL1L1I)
     YLLI = TEMP1*(COSPH*YL1L1I + SINPH*YL1L1R)
     Y(INDEX2) = DCMPLX(YLLR,YLLI)
     Y(INDEX)  = MSIGN*DCONJG(Y(INDEX2))
     INDEX2 = INDEX2 - 1
     INDEX  = INDEX  + 1

!    YLL1 = Y(L,L-1) = f(Y(L-1,L-1)) ... Formula 2
!           (the coefficient for Y(L-2,L-1) in Formula 2 is zero)
     TEMP2 = sqrt(dble(2*L+1))*COSTH
     YLL1R = TEMP2*YL1L1R
     YLL1I = TEMP2*YL1L1I
     Y(INDEX2) = DCMPLX(YLL1R,YLL1I)
     Y(INDEX)  = -MSIGN*DCONJG(Y(INDEX2))
     INDEX2 = INDEX2 - 1
     INDEX  = INDEX  + 1

     I4L2 = INDEX2 - 4*L + 2
     I2L  = INDEX2 - 2*L
     D4LL1C = COSTH*sqrt(dble(4*L*L-1))
     D2L13  = -sqrt(dble(2*L+1)/dble(2*L-3))

     do M = L - 2, 0, -1
!       YLM = Y(L,M) = f(Y(L-2,M),Y(L-1,M)) ... Formula 2
        TEMP1 = 1.0D+0/sqrt(dble((L+M)*(L-M)))
        TEMP2 = D4LL1C*TEMP1
        TEMP3 = D2L13*sqrt(dble((L+M-1)*(L-M-1)))*TEMP1
        YLMR = TEMP2*dble(Y(I2L))  + TEMP3*dble(Y(I4L2))
        YLMI = TEMP2*aimag(Y(I2L)) + TEMP3*aimag(Y(I4L2))
        Y(INDEX2) = DCMPLX(YLMR,YLMI)
        Y(INDEX)  = MSIGN*DCONJG(Y(INDEX2))

        MSIGN  = -MSIGN
        INDEX2 = INDEX2 - 1
        INDEX  = INDEX  + 1
        I4L2   = I4L2   - 1
        I2L    = I2L    - 1
     end do
  end do
end subroutine Ylm
end module     Ylm_m

!!/---
!! Local Variables:
!! mode: f90
!! End:
!!\---
!!
!! Time-stamp: <2016-08-04 09:44:17 assman@faepop71.tu-graz.ac.at>
