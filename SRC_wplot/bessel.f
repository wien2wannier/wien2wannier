!!! wien2wannier/SRC_wplot/bessel.f

SUBROUTINE BESSEL(NK,LDNK,BK,NMT,RAD,LMAX2,F,DF)
  use const
  IMPLICIT REAL(R8) (A-H,O-Z)
  DIMENSION BK(3,NK), RAD(NMT)
  DIMENSION F(LMAX2+1,LDNK,NMT), DF(LMAX2+1,LDNK,NMT)
  !
  LDIM = LMAX2 + 1
  DO IR=1,NMT
     RMT = RAD(IR)
     DO IK=1,NK
        AK  = SQRT( BK(1,IK)*BK(1,IK) + BK(2,IK)*BK(2,IK) + &
             BK(3,IK)*BK(3,IK) )
        ARG = AK * RMT
        CALL SPHBES(LMAX2,ARG,F(1,IK,IR))
        CALL DVBES1(F(1,IK,IR),DF(1,IK,IR),ARG,RMT,LDIM)
        DO L=1,LDIM
           DF(L,IK,IR) = DF(L,IK,IR)*AK
        END DO
     END DO
  END DO
END SUBROUTINE BESSEL


!!/---
!! Local Variables:
!! mode: f90
!! End:
!!\---
!!
!! Time-stamp: <2015-05-23 19:58:48 elias>
