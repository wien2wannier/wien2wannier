!!! wien2wannier/SRC_wplot/sphbes.f

subroutine SPHBES(N,X,FJ)
  use const, only: DPk

  implicit none

  integer,   intent(in)  :: N
  real(DPk), intent(in)  :: X
  real(DPk), intent(out) :: FJ(*)

  real(DPk), parameter :: &
       XLIM=0.1_DPk, HF=0.5_DPk, TNHF=10.5_DPk, FFT=15.0_DPk, &
       T25=1e25_DPk, TN25=1e-25_DPk, TN50=1e-50_DPk

  real(DPk) :: HFXSQ, XL, TWM, TA, XLP, CUFAC, FFO, FFN, XI, FM, FFP
  real(DPk) :: SDR, SER
  integer   :: m, mm, j, jj, NS

!***********************************************************************
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
!***********************************************************************
      IF(N < 0) STOP 'SPHBES: N SHOULD NOT BE NEGATIVE'
      IF(X < 0) STOP 'SPHBES: X SHOULD NOT BE NEGATIVE'
      IF (X > XLIM) GO TO 25
      HFXSQ=HF*X*X
      XL=1
      TWM=1
      M=0
  11  M=M+1
      TA=XL
      TWM=TWM+2
      XL=XL/TWM
      TA=TA-XL*HFXSQ
      XLP=XL/(TWM+2)
      FJ(M)=TA+HF*XLP*HFXSQ*HFXSQ
      XL=XL*X
      IF (M.LE.N)  GO TO 11
      RETURN
  25  CUFAC=4.2_DPk
      IF (X.LT.(N-2)) CUFAC=TNHF/(N+HF-X)
      NS = N+5 + int(X*CUFAC)
!*******************  ADD ADDITIONAL FACTOR  ***************************
      NS=NS + int(FFT/(1+SQRT(X)))
!***********************************************************************
      CONTINUE
      FFO=0
      FFN=TN25
      M=NS-1
      XI=1/X
      FM=(M+M)+1
      SDR=FM*TN50
 314  FFP=FM*XI*FFN-FFO
      IF (ABS(FFP).LT.T25) GO TO 315
      SDR=SDR*TN50
      FFP=FFP*TN25
      FFN=FFN*TN25
 315  SDR=SDR + (FM-2)*FFP*FFP
      FFO=FFN
      FFN=FFP
      IF (M.LE.N) GO TO 316
      M=M-1
      FM=FM-2
      GO TO 314
 316  FJ(M)=FFN
      FJ(M+1)=FFO
      GO TO 33
  32  FJ(M)=FM*XI*FJ(M+1)-FJ(M+2)
      IF(ABS(FJ(M)).GE.T25) GO TO 56
      SDR=SDR + (FM-2)*FJ(M)*FJ(M)
      IF (M.LE.1) GO TO 34
   33 M = M-1
      FM=FM-2
      GO TO 32
 34   SER=1/SQRT(SDR)
      MM = N+1
      DO 40 M=1,MM
      FJ(M)=FJ(M)*SER
  40  CONTINUE
      RETURN
   56 JJ= M+1
      NS=N+1
      DO 57 J = JJ,NS
      FJ(J)=FJ(J)*TN25
  57  CONTINUE
      SDR=SDR*TN50
      GO TO 32
      END


!!/---
!! Local Variables:
!! mode: f90
!! End:
!!\---
!!
!! Time-stamp: <2016-07-15 11:46:00 assman@faepop71.tu-graz.ac.at>
