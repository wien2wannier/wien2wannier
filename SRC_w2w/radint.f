!!! wien2wannier/SRC_w2w/radint.f

      SUBROUTINE RADINT(JATOM,LJMAX,BM)
      USE param
      USE struct
      USE bessel
      use lolog, only: n_rad
      use const, only: R8
      use radfu, only: RF1, RF2
      use gener, only: br1

      IMPLICIT REAL(R8) (A-H,O-Z)
      INTEGER L_index,R_index

      DIMENSION  A(nrad),B(nrad),X(nrad),Y(nrad)


      DO  I=1,JRJ(JATOM)
       RX=R0(JATOM)*EXP(DX(JATOM)*(i-1))*BM
       call sphbes(LJMAX+1,RX,rj(0,I))
      ENDDO
      L_index=0
      DO 19 L1=0,LMAX2
         DO 20 L2=0,LMAX2
            DO 21 LJ=0,LJMAX
       IF (MOD((L1+L2+LJ),2) .EQ. 1) GOTO 280
       IF ((L1+L2-LJ) .LT. 0) GOTO 280
       IF ((L1-L2+LJ) .LT. 0) GOTO 280
       IF ((-L1+L2+LJ) .LT. 0) GOTO 280
       L_index=L_index+1
       R_index=0
       DO 30 IF1=1,n_rad(l1)
       DO 31 IF2=1,n_rad(l2)
        R_index=R_index+1
        DO  I=1,JRJ(JATOM)
         A(i)=rf1(i,l1,if1)*rj(lj,i)
         B(i)=rf2(i,l1,if1)*rj(lj,i)
         X(i)=rf1(i,l2,if2)
         Y(i)=rf2(i,l2,if2)
        ENDDO
        CALL RINT13(A,B,X,Y, &
        ri_mat(r_index,l_index),JATOM)

 31   CONTINUE
 30   CONTINUE
 280  CONTINUE
 21   CONTINUE
 20   CONTINUE
 19   CONTINUE
      END


!!/---
!! Local Variables:
!! mode: f90
!! End:
!!\---
!!
!! Time-stamp: <2016-07-05 15:09:34 assman@faepop71.tu-graz.ac.at>
