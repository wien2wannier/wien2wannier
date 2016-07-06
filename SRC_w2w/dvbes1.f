!!! wien2wannier/SRC_w2w/dvbes1.f

      SUBROUTINE DVBES1(FJ,DJ,SM,NT)
!-----X----X----X----X----X----X----X----X----X----X----X----X----X----X
!-----X CALCULATE THE DERIVATIVES OF THE BESSEL FUNCTIONS.   X----X----X
!-----X   DJ=DFJ/DX WHERE X=SM*RI                                 X----X
!-----X                    D.D.KOELLING                      X----X----X
!-----X----X----X----X----X----X----X----X----X----X----X----X----X----X
      use const, only: R8
      implicit none

      real(R8) :: FJ(NT), DJ(NT), SM
      integer  :: NT

      intent(in)  :: FJ, SM, NT
      intent(out) :: DJ

      real(R8), parameter :: ZUP = 1.0e-5_R8

      integer  :: l, lm
      real(R8) :: q2, q3

      IF(SM.GT.ZUP) GOTO 20
      DJ(1)=0
      DJ(2)=1/3._R8
      DO 10 L=3,NT
         DJ(L)=0
 10   CONTINUE
      RETURN
   20 Q2=-1/SM
      Q3=Q2
      DJ(1)=-FJ(2)
      LM=1
      DO 30 L=2,NT
      Q3=Q3+Q2
      DJ(L)=FJ(LM)+Q3*FJ(L)
      LM=LM+1
   30 CONTINUE
      RETURN
      END


!!/---
!! Local Variables:
!! mode: f90
!! End:
!!\---
!!
!! Time-stamp: <2016-07-05 16:15:19 assman@faepop71.tu-graz.ac.at>
