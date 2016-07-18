!!! wien2wannier/SRC_w2w/rotate.f

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

!!/---
!! Local Variables:
!! mode: f90
!! End:
!!\---
!!
!! Time-stamp: <2016-07-18 15:50:58 assman@faepop71.tu-graz.ac.at>
