!!! wien2wannier/SRC_w2w/harmon.f

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
  DO I=1,N
     A(1)=X(I)*BR1(1,1)+Y(I)*BR1(1,2)+Z(I)*BR1(1,3)
     A(2)=X(I)*BR1(2,1)+Y(I)*BR1(2,2)+Z(I)*BR1(2,3)
     A(3)=X(I)*BR1(3,1)+Y(I)*BR1(3,2)+Z(I)*BR1(3,3)
     XM=SQRT(A(1)**2+A(2)**2+A(3)**2)
     XA=RI*XM
     CALL SPHBES(LMAX2,XA,F(1,I))
     CALL DVBES1(F(1,I),DF(1,I),XA,LMX)
     do J=1,LMX
        DF(J,I)=XM*DF(J,I)
     end do
  end do
end subroutine harmon
end module harmon_m

!!/---
!! Local Variables:
!! mode: f90
!! End:
!!\---
!!
!! Time-stamp: <2016-07-18 15:47:34 assman@faepop71.tu-graz.ac.at>
