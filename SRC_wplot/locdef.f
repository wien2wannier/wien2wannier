!!! wien2wannier/SRC_wplot/locdef.f

SUBROUTINE LOCDEF(ROT0,IMAT,ROT)
  use const
  use latt

  !:17[
  use param
  !:17]

  IMPLICIT  REAL(R8) (A-H,O-Z)
  DIMENSION ROT0(3,3),IMAT(3,3),ROT(3,3)
  !
  ! Input:
  ! ROT0 : local reference rotation matrix (Cartesian coordinates)
  !        x'_i = Sum(j) (T^-1)_ij x_j  with  T_ji = (T^-1)_ij = ROT0(i,j)
  ! IMAT : symmetry operation {Q|t} (primitive fractional coordinates)
  !        y_m = Sum(n) Q_mn x_n + t_m  with  Q_mn = IMAT(n,m)
  !
  ! from module LATT
  ! BR2  : primitive real space lattice vector am = BR2(m,:)
  ! BR4  : primitive reciprocal lattice vector bn = BR4(n,:) [without 2pi]
  !
  ! Output:
  ! ROT  : symmetry adapted local rotation matrix (Cartesian coordinates)
  !        x'_i = Sum(j) (R^-1)_ij x_j  with  R_ji = (R^-1)_ij = ROT(i,j)
  !        [ ROT and ROT0 may share memory ! ]
  !
  ! Algorithm:
  ! R := Q o T  with  Q_ik = Sum(m,n) am_i Q_mn bn_k

  DIMENSION AMAT(3,3),BMAT(3,3)
  !
  !:17[
  IF(.NOT.USEROT)THEN
     !       << don't apply any local rotations, i.e. R := E >>
     ROT = 0

     DO K=1,3
        ROT(K,K)=1.0D0
     END DO
     RETURN
  ENDIF
  !:17]
  !     << A_mk := Sum(n) Q_mn bn_k >>
  DO M=1,3
     DO K=1,3
        AMAT(M,K)=IMAT(1,M)*BR4(1,K) &
             +IMAT(2,M)*BR4(2,K) &
             +IMAT(3,M)*BR4(3,K)
     END DO
  END DO
  !
  !:17[
  IF(.NOT.ADDLOC)THEN
     !       << ignore local rotation matrix T from input, i.e. R = Q >>
     !       << R_ik = Q_ik = Sum(m) am_i A_mk >>
     DO I=1,3
        DO K=1,3
           ROT(K,I)=BR2(1,I)*AMAT(1,K) &
                +BR2(2,I)*AMAT(2,K) &
                +BR2(3,I)*AMAT(3,K)
        END DO
     END DO
     RETURN
  ENDIF
  !:17]
  !     << B_mj = Sum(k) A_mk T_kj >>
  DO M=1,3
     DO J=1,3
        BMAT(M,J)=AMAT(M,1)*ROT0(J,1) &
             +AMAT(M,2)*ROT0(J,2) &
             +AMAT(M,3)*ROT0(J,3)
     END DO
  END DO
  !
  !     << R_ij = Sum(m) am_i B_mj >>
  DO I=1,3
     DO J=1,3
        ROT(J,I)=BR2(1,I)*BMAT(1,J) &
             +BR2(2,I)*BMAT(2,J) &
             +BR2(3,I)*BMAT(3,J)
     END DO
  END DO
  !
  RETURN
END SUBROUTINE LOCDEF


!!/---
!! Local Variables:
!! mode: f90
!! End:
!!\---
!!
!! Time-stamp: <2015-05-23 19:58:48 elias>
