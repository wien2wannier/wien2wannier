SUBROUTINE RADINT(JATOM,LJMAX,BM)
USE param
USE struct
USE bessel
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER L_index,R_index
      COMMON /RADFU/   RF1(NRAD,0:LMAX2,nrf),RF2(NRAD,0:LMAX2,nrf)
      COMMON /GENER/  BR1(3,3)
      logical loor(0:lomax),lapw(0:lmax2)
      common /lolog/   nlo,nlov,nlon,loor,ilo(0:lomax),lapw,n_rad(0:lmax2)
      DIMENSION  A(nrad),B(nrad),X(nrad),Y(nrad)


      DO  I=1,JRJ(JATOM)
       RX=R0(JATOM)*EXP(DX(JATOM)*(i-1))*BM
       call sphbes(LJMAX,RX,rj(0,I))
!       write(75,5)R0(JATOM)*EXP(DX(JATOM)*(i-1)),rf1(I,0,1),rf1(I,0,2),rf1(I,0,3)
 5     format(10f12.8)
      ENDDO
      L_index=0
      DO 19 L1=0,LMAX2
      DO 19 L2=0,LMAX2
      DO 19 LJ=0,LJMAX
       IF (MOD((L1+L2+LJ),2) .EQ. 1) GOTO 280
       IF ((L1+L2-LJ) .LT. 0) GOTO 280
       IF ((L1-L2+LJ) .LT. 0) GOTO 280
       IF ((-L1+L2+LJ) .LT. 0) GOTO 280
       L_index=L_index+1
       R_index=0
       DO 20 IF1=1,n_rad(l1)
       DO 20 IF2=1,n_rad(l2)
        R_index=R_index+1
        DO  I=1,JRJ(JATOM)
         A(i)=rf1(i,l1,if1)*rj(lj,i)
         B(i)=rf2(i,l1,if1)*rj(lj,i)
         X(i)=rf1(i,l2,if2)
         Y(i)=rf2(i,l2,if2)
        ENDDO
        CALL RINT13(A,B,X,Y, &
        ri_mat(r_index,l_index),JATOM)
!       write(76,6)l_index,l2,lj,l1,if2,if1,ri_mat(r_index,l_index)
  6     format(i3,4x,5i3,f14.9)
 20    CONTINUE
 280   CONTINUE
 19   CONTINUE
      END

