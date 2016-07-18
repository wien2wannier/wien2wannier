!!! wien2wannier/SRC_w2w/rint13.f

module     rint13_m; contains
subroutine rint13(A,B,X,Y,S,JATOM)
  !     PERFORM RADIAL INTEGRALS REQUIRED BY BHDK13
  !                            D.D.KOELLING

  use param,  only: clight, Nrad
  use struct, only: rel, dx, jrj, R0
  use const,  only: R8

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

!!/---
!! Local Variables:
!! mode: f90
!! End:
!!\---
!!
!! Time-stamp: <2016-07-18 15:51:07 assman@faepop71.tu-graz.ac.at>
