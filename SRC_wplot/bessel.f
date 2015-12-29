!!! wien2wannier/SRC_wplot/bessel.f

SUBROUTINE BESSEL(NK,LDNK,BK,NMT,RAD,LMAX2,F,DF)
  use param, only: DPk

  implicit none

  integer,   intent(in)  :: Nk, LDNK, lMax2, NMT
  real(DPk), intent(in)  :: BK(3,Nk), rad(NMT)
  real(DPk), intent(out) :: F(LMAX2+1,LDNK,NMT), DF(LMAX2+1,LDNK,NMT)

  integer   :: ldim, ir, ik, l
  real(DPk) :: AK, arg, RMT

  LDIM = LMAX2 + 1
  DO IR=1,NMT
     RMT = RAD(IR)
     DO IK=1,NK
        AK  = SQRT( BK(1,IK)*BK(1,IK) + BK(2,IK)*BK(2,IK) + &
             BK(3,IK)*BK(3,IK) )
        ARG = AK * RMT
        CALL SPHBES(LMAX2,ARG,F(1,IK,IR))
        CALL DVBES1(F(:,IK,IR),DF(:,IK,IR),ARG,LDIM)
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
!! Time-stamp: <2015-12-29 16:25:52 assman@faepop36.tu-graz.ac.at>
