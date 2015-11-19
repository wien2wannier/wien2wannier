      SUBROUTINE BESSEL(NK,LDNK,BK,NMT,RAD,LMAX2,F,DF)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION BK(3,NK), RAD(NMT)
      DIMENSION F(LMAX2+1,LDNK,NMT), DF(LMAX2+1,LDNK,NMT)
!
      LDIM = LMAX2 + 1
      DO 10 IR=1,NMT
        RMT = RAD(IR)
        DO 20 IK=1,NK
          AK  = SQRT( BK(1,IK)*BK(1,IK) + BK(2,IK)*BK(2,IK) + &
                      BK(3,IK)*BK(3,IK) )
          ARG = AK * RMT
          CALL SPHBES(LMAX2,ARG,F(1,IK,IR))
          CALL DVBES1(F(1,IK,IR),DF(1,IK,IR),ARG,RMT,LDIM)
          DO 30 L=1,LDIM
            DF(L,IK,IR) = DF(L,IK,IR)*AK
   30     CONTINUE
   20   CONTINUE
   10 CONTINUE
!
      RETURN
      END
