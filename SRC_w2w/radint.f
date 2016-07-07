!!! wien2wannier/SRC_w2w/radint.f

subroutine RADINT(JATOM,LJMAX,BM)
  use param,  only: Lmax2, Nrad
  use struct, only: jrj, R0, dx
  use bessel, only: rj, ri_mat
  use lolog,  only: n_rad
  use const,  only: R8
  use radfu,  only: RF1, RF2

  implicit none

  integer,  intent(in) :: jatom, LJmax
  real(R8), intent(in) :: bm

  real(R8) :: A(Nrad), B(Nrad), X(Nrad), Y(Nrad), RX
  integer  :: L_index,l1,l2,lj, R_index, i, if1,if2

  DO  I=1,JRJ(JATOM)
     RX=R0(JATOM)*EXP(DX(JATOM)*(i-1))*BM
     call sphbes(LJMAX+1,RX,rj(0,I))
  ENDDO
  L_index=0
  do L1=0,LMAX2
     DO L2=0,LMAX2
        ljloop: do LJ=0,LJMAX
           if (mod((L1+L2+LJ),2) == 1  &
                .or. (L1+L2-LJ)   < 0  &
                .or. (L1-L2+LJ)   < 0  &
                .or. (-L1+L2+LJ)  < 0) &
                cycle ljloop
           L_index=L_index+1
           R_index=0
           do IF1=1,n_rad(l1)
              do IF2=1,n_rad(l2)
                 R_index=R_index+1
                 do  I=1,JRJ(JATOM)
                    A(i)=rf1(i,l1,if1)*rj(lj,i)
                    B(i)=rf2(i,l1,if1)*rj(lj,i)
                    X(i)=rf1(i,l2,if2)
                    Y(i)=rf2(i,l2,if2)
                 enddo
                 call RINT13(A,B,X,Y, &
                      ri_mat(r_index,l_index),JATOM)

              end do
           end do
        end do ljloop
     end do
  end do
end subroutine RADINT


!!/---
!! Local Variables:
!! mode: f90
!! End:
!!\---
!!
!! Time-stamp: <2016-07-07 10:18:18 assman@faepop71.tu-graz.ac.at>
