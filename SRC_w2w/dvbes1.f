!!! wien2wannier/SRC_w2w/dvbes1.f

subroutine DVBES1(FJ,DJ,SM,NT)
!-----X----X----X----X----X----X----X----X----X----X----X----X----X----X
!-----X CALCULATE THE DERIVATIVES OF THE BESSEL FUNCTIONS.   X----X----X
!-----X   DJ=DFJ/DX WHERE X=SM*RI                                 X----X
!-----X                    D.D.KOELLING                      X----X----X
!-----X----X----X----X----X----X----X----X----X----X----X----X----X----X
  use const, only: DPk

  implicit none

  real(DPk), intent(in)  :: FJ(*), SM
  real(DPk), intent(out) :: DJ(*)
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
end subroutine DVBES1


!!/---
!! Local Variables:
!! mode: f90
!! End:
!!\---
!!
!! Time-stamp: <2016-07-15 12:09:17 assman@faepop71.tu-graz.ac.at>
