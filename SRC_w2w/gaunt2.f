!!! wien2wannier/SRC_w2w/gaunt2.f

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
     STH = SQRT(1-CTH*CTH)
     RF = 1/SQRT(FPI)
     YR(K,1) = RF*FACTOR
     I = 1
     P(1,1) = 1
     P(2,1) = CTH
     C2L = CTH
     TCTH = CTH + CTH
     L1 = 2
     L = 1
10   CONTINUE
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
20   CONTINUE
     !
     !        recurse upward in L
     !
     P(L1,M) = (C2L*P(L,M)-LM*P(L2,M))/LM2
     C1L = (LM+1)*CTH
     P(L,M1) = 0
     !
     !        recurse upward in M
     !
     IF (ABS(STH) .GE. SNULL) P(L,M1) = (C1L*P(L,M)-LM2*P(L1,M))/STH
30   CONTINUE
     I = I + 1
     IDWN = IDWN - 1
     CSR = SQRT((2*L-1)/(FPI*CD))
     CYP = SGNM*CSR*P(L,M)
     YR(K,I) = CYP*FACTOR
     IX = I - (L*L - L + 1)
     IF (IDWN .NE. I) YR(K,IDWN) = FACTOR*CYP*(-1)**IX
     M = M1
     M1 = M + 1
     LM2 = LM2 - 1
     LM = LM + 1
     CD = CD*LM*LM2
     SGNM = -SGNM
     IF (M - L < 0) goto 20
     IF (M - L ==0) goto 30
     IF (L .LE. LOMAX) GOTO 10
  end do
end subroutine gaunt2
end module gaunt2_m


!!/---
!! Local Variables:
!! mode: f90
!! End:
!!\---
!!
!! Time-stamp: <2016-07-18 13:57:13 assman@faepop71.tu-graz.ac.at>
