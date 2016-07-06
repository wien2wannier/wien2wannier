!!! wien2wannier/SRC_w2w/dergl.f

!rschmid
!
!     Calculate the first derivate of f_l.
!
!rschmid

subroutine dergl(g_l,der1,rnot,dx,mesh)
  use param
  use const, only: R8

  implicit none

  real(R8) :: g_l(nrad), der1(nrad), dx, rnot
  integer  :: mesh

  intent(in)  :: g_l, dx, rnot, mesh
  intent(out) :: der1

  real(R8) :: f_l(nrad), dfldi(nrad), r(nrad), rad
  integer  :: j,i

  do i=1,mesh
     r(i) = rnot*exp(dx*(i-1))
     f_l(i) = g_l(i)/r(i)
  enddo

  !rschmid
  !  Use unsymmetric 6-point formulae for the first derivative at the
  !  first 3 mesh points.
  !rschmid
  dfldi(1) =  ( - 137.D0*f_l(1) + 300.D0*f_l(2)    &
       - 300.D0*f_l(3) + 200.D0*f_l(4) &
       -  75.D0*f_l(5) +  12.D0*f_l(6) ) &
       / 60.D0
  dfldi(2) = ( -  12.D0*f_l(1) -  65.D0*f_l(2) &
       + 120.D0*f_l(3) -  60.D0*f_l(4) &
       +  20.D0*f_l(5) -   3.D0*f_l(6) ) &
       / 60.D0
  dfldi(3) = (     3.D0*f_l(1) -  30.D0*f_l(2) &
       -  20.D0*f_l(3) +  60.D0*f_l(4)  &
       -  15.D0*f_l(5) +   2.D0*f_l(6) ) &
       / 60.D0
  RAD=R(1)
  DER1(1)=dfldi(1)/RAD/DX
  RAD=R(2)
  DER1(2)=dfldi(2)/RAD/DX
  RAD=R(3)
  DER1(3)=dfldi(3)/RAD/DX
  !rschmid
  !  Use symmetric 7-point formula to generate the first derivative at
  !  all intermediate mesh points.
  !rschmid

  do J=4,MESH-3
     dfldi(J) = (        f_l(J+3) - f_l(J-3) &
          -  9.D0 * (f_l(J+2) - f_l(J-2))  &
          + 45.D0 * (f_l(J+1) - f_l(J-1)) ) &
          / 60.D0
     RAD=R(J)
     DER1(J)=dfldi(J)/RAD/DX
  end do
!rschmid
!  Use unsymmetric 6-point formulae for the second derivative at the
!  last 3 mesh points.
!rschmid
  dfldi(MESH-2) = - (    3.D0*f_l(MESH  ) &
       -  30.D0*f_l(MESH-1) &
       -  20.D0*f_l(MESH-2) &
       +  60.D0*f_l(MESH-3)   &
       -  15.D0*f_l(MESH-4)   &
       +   2.D0*f_l(MESH-5) ) &
       / 60.D0
  RAD=R(MESH-2)
  DER1(MESH-2)=dfldi(MESH-2)/RAD/DX

  dfldi(MESH-1) = - (-  12.D0*f_l(MESH  ) &
       -  65.D0*f_l(MESH-1) &
       + 120.D0*f_l(MESH-2) &
       -  60.D0*f_l(MESH-3)   &
       +  20.D0*f_l(MESH-4)   &
       -   3.D0*f_l(MESH-5) ) &
       / 60.D0
  RAD=R(MESH-1)
  DER1(MESH-1)=dfldi(MESH-1)/RAD/DX

  dfldi(MESH  ) = - (- 137.D0*f_l(MESH  ) &
       + 300.D0*f_l(MESH-1) &
       - 300.D0*f_l(MESH-2) &
       + 200.D0*f_l(MESH-3) &
       -  75.D0*f_l(MESH-4) &
       +  12.D0*f_l(MESH-5) ) &
       / 60.D0

  RAD=R(MESH)
  DER1(MESH)=dfldi(MESH)/RAD/DX
end subroutine dergl

!!/---
!! Local Variables:
!! mode: f90
!! End:
!!\---
!!
!! Time-stamp: <2016-07-06 16:33:28 assman@faepop71.tu-graz.ac.at>
