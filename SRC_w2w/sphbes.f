!!! wien2wannier/SRC_w2w/sphbes.f

subroutine SPHBES(N,X,FJ)
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

  use const, only: R8
  use clio,  only: croak
  use util,  only: string

  implicit none

  integer,  intent(in)  :: N
  real(R8), intent(in)  :: X
  real(R8), intent(out) :: FJ(*)

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
end subroutine SPHBES


!!/---
!! Local Variables:
!! mode: f90
!! End:
!!\---
!!
!! Time-stamp: <2016-07-15 15:12:39 assman@faepop71.tu-graz.ac.at>
