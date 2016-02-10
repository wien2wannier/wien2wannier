!!! wien2wannier/SRC_wplot/trans.f

subroutine trans(stru)
  use param,     only: DPk
  use latt,      only: br1, br2, br3, br4
  use param,     only: Nsym
  use sym2,      only: tau, imat
  use structmod, only: struct_t

  implicit none

  type(struct_t), intent(inout) :: stru

  real(DPk) :: F(3), S(3,3), T(3,3), Q(3,3)
  integer   :: i, k, n
!
! transforms real space vectors x and symmetry operations {Q|t}
! from conventional into primitive fractional coordinates
! --------------------------------------------------------------------
! Input:
! NPOS     -- the number of real space vectors x to transform
! NSYM     -- the number of symmetry operations {Q|t} to transform
!
! module LATT
! BR1(i,:) -- conventional real space lattice vectors a_i
! BR2(i,:) -- primitive    real space lattice vectors a_i
! BR3(i,:) -- conventional reciprocal lattice vectors b_i (without 2pi)
! BR4(i,:) -- primitive    reciprocal lattice vectors b_i (without 2pi)
!
! Input/Output:
! POS(:,n)    -- the n-th real space vector x
!                x = Sum(i=1,3) x_i a_i  with  x_i = POS(i,n)
! IMAT(:,:,n) -- the n-th symmetry operation {Q,t} :
! TAU(:,n)       y_i = Sum(j) Q_ij x_i + t_i  with  Q_ij = IMAT(j,i,n)
!                                             and   t_i  = TAUK(  i,n)
!
! Algorithm:
! real space vectors from conventional to primitive:
! p_k = Sum(i) T(k,i) c_i  with  T(k,i) = Sum(j) BR4(k,j) BR1(i,j)
!
! real space vectors from primitive to conventional:
! c_i = Sum(k) S(i,k) p_k  with  S(i,k) = Sum(j) BR3(i,j) BR2(k,j)
!
! symmetry operations from conventional to primitive
! Q_kk' = Sum(ii') T(k,i) Q_ii' S(i',k')  and  t_k = Sum(i) T(k,i) t_i
! --------------------------------------------------------------------
!

!     << set up the transformation matrices >>
  DO K=1,3
     DO I=1,3
        T(K,I) = BR4(K,1)*BR1(I,1) + BR4(K,2)*BR1(I,2) &
             + BR4(K,3)*BR1(I,3)
        S(I,K) = BR3(I,1)*BR2(K,1) + BR3(I,2)*BR2(K,2) &
             + BR3(I,3)*BR2(K,3)
     END DO
  END DO

!     << transform the real space vectors >>
  DO N=1,stru%Nat
     DO K=1,3
        F(K) = T(K,1)*stru%pos(1,N) + T(K,2)*stru%pos(2,N) &
             + T(K,3)*stru%pos(3,N)
     END DO
     DO K=1,3
        stru%pos(K,N) = F(K)
     END DO
  END DO

!     << transform the symmetry operations >>
  DO N=1,NSYM
     DO K=1,3
        F(K) = T(K,1)*TAU(1,N) + T(K,2)*TAU(2,N) &
             + T(K,3)*TAU(3,N)
        DO I=1,3
           Q(K,I) = T(K,1)*IMAT(I,1,N) + T(K,2)*IMAT(I,2,N) &
                & + T(K,3)*IMAT(I,3,N)
        END DO
     END DO
     DO K=1,3
        TAU(K,N) = F(K)
        DO I=1,3
           IMAT(I,K,N) = NINT( Q(K,1)*S(1,I) + Q(K,2)*S(2,I) + &
                &              Q(K,3)*S(3,I) )
        END DO
     END DO
  END DO

END SUBROUTINE TRANS


!!/---
!! Local Variables:
!! mode: f90
!! End:
!!\---
!!
!! Time-stamp: <2016-02-10 16:59:04 assman@faepop36.tu-graz.ac.at>
