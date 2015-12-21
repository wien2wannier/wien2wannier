!!! wien2wannier/SRC_w2w/planew.F

      SUBROUTINE PLANEW(NB,num_kpts,NNTOT,maxx,maxy,maxz)     
  USE param
  USE struct
  USE xa
  USE xa3
  USE amn_mmn
  USE pairs
      IMPLICIT REAL(R8) (A-H,O-Z)
      logical          loor(0:lomax),lapw(0:lomax)
      logical, allocatable, dimension(:,:) :: a_null,a_null_fin
      INTEGER          pair_index
      COMPLEX(C16)       SUM(NB),SUMA(NB,NB),FF
      complex(c16), allocatable, dimension(:,:,:) :: FAC
      COMPLEX(C16)       PHSHEL
      COMPLEX(C16), allocatable, dimension(:,:,:) :: A_
      COMPLEX(C16), allocatable, dimension(:,:) :: B_
 
      COMMON /GENER/  BR1(3,3),BR2(3,3)                                 
      common /lolog/  nlo,nlov,nlon,loor,ilo(0:lomax),lapw  

      !allocate memory to arrays
      allocate(FAC(-(2*maxx+1):(2*maxx+1),-(2*maxy+1):(2*maxy+1),-(2*maxz+1):(2*maxz+1)) )
      allocate(A_(NB,0:(2*maxx+1)*(2*maxy+1)*(2*maxz+1)-1,NUM_KPTS) )
      allocate(B_(NB,0:(2*maxx+1)*(2*maxy+1)*(2*maxz+1)-1) )
      allocate(a_null(0:(2*maxx+1)*(2*maxy+1)*(2*maxz+1)-1,NUM_KPTS) )     
      allocate(a_null_fin(0:(2*maxx+1)*(2*maxy+1)*(2*maxz+1)-1,NUM_KPTS) )
    
      TWOPI=2.D0*PI               
      m2x=2*maxx+1     
      m2y=2*maxy+1  
      m2z=2*maxz+1 
      max=m2x*m2y*m2z-1 
      write(unit_out,*)'MAX Q:',maxx,maxy,maxz
     t_b=0.d0
     t_prod=0.d0
       
      call cputim(t1)  
      A_=0
      a_null=.false.
      a_null_fin=.false.
      pair_index=0
      DO k1=1,num_kpts
       pair_index=pair_index+1
        kkk=KP(pair_index)
        DO I=1,size(kkk)-(nlo+nlon+nlov)
        I_X=GX(I,kkk)+maxx
        I_Y=GY(I,kkk)+maxy
        I_Z=GZ(I,kkk)+maxz
        index=(I_Z*m2y+I_Y)*m2x+I_X
        a_null(index,kkk)=.true.
        DO num=1,nb
           A_(num,index,kkk)=A(I,num,kkk)
        ENDDO
       ENDDO
       DO k2=1,NNTOT
        if (k2.gt.1) pair_index=pair_index+1
        kk=KPB(pair_index)
        DO I=1,size(kkk)
         I_X=GX(I,kk)+maxx-BQX(pair_index)
         I_Y=GY(I,kk)+maxy-BQY(pair_index)
         I_Z=GZ(I,kk)+maxz-BQZ(pair_index)
         index=(I_Z*m2y+I_Y)*m2x+I_X
         a_null_fin(index,kkk)=.true.
        ENDDO
       ENDDO
      ENDDO       
      
      FAC=0
      LFIRST=1
      DO JATOM=1,NAT
       IF (JATOM.GT.1) LFIRST=LFIRST+MULT(JATOM-1)
       VV=4*PI*RMT(JATOM)**3/3
       DO K=-m2z,m2z
        DO J=-m2y,m2y
         DO I=-m2x,m2x
          BX=I*BR1(1,1)+J*BR1(1,2)+K*BR1(1,3)
          BY=I*BR1(2,1)+J*BR1(2,2)+K*BR1(2,3)
          BZ=I*BR1(3,1)+J*BR1(3,2)+K*BR1(3,3)
          X=RMT(JATOM)*SQRT(BX*BX+BY*BY+BZ*BZ)
          PHSHEL=0
          DO MU=1,MULT(JATOM)
           LATOM=LFIRST+MU-1
           ARG=I*POS(1,LATOM)  &
              +J*POS(2,LATOM)  &
              +K*POS(3,LATOM)
           PHSHEL=PHSHEL+exp(-(0,1)*ARG*TWOPI)
          ENDDO
          IF ((I*I+J*J+K*K).eq.0) THEN
           FAC(I,J,K)=FAC(I,J,K)-MULT(JATOM)*VV/VOL
          ELSE
           FAC(I,J,K)=FAC(I,J,K)-3*PHSHEL*(SIN(X)-X*COS(X))/(X*X*X)*VV/VOL
          ENDIF
         ENDDO
        ENDDO
       ENDDO
      ENDDO
      FAC(0,0,0)=FAC(0,0,0)+1.d0
      call cputim(t2)
      write(unit_out,*)'t_init=',t2-t1

!............................................................................
!.....READ IN K POINT AND BASIS VECTORS  
      pair_index=0
      DO 3 k1=1,num_kpts
       pair_index=pair_index+1
       kkk=KP(pair_index)
       call cputim(tt1)
       b_=0
       DO J=0,max
        if (a_null_fin(J,kkk)) then
         SUM=0
         JJ=J/m2x
         iz2=jj/m2y
         iy2=jj-iz2*m2y
         ix2=j-jj*m2x
         DO I=0,max
          if(a_null(I,kkk)) then
           II=I/m2x
           iz1=II/m2y
           iy1=II-iz1*m2y
           ix1=I-II*m2x
           FF=FAC(ix2-ix1,iy2-iy1,iz2-iz1)
           sum=sum+a_(1:nb,I,k1)*FF
          endif
         ENDDO
         b_(1:nb,j)=sum
        endif
       ENDDO
       call cputim(tt2)
       DO 4 k2=1,NNTOT
        if (k2.gt.1) pair_index=pair_index+1
        kkk=KPB(pair_index)
        !bugfix 07/23 p.wissgott 
   !    id=-(BQZ(pair_index)*m2y+BQY(pair_index))*m2z-BQX(pair_index)-1
        id=-(BQZ(pair_index)*m2y+BQY(pair_index))*m2x-BQX(pair_index)-1
        SUMA=0
!debug
 !id=0
        DO J=0,max   
         id=id+1
!        id=j-(BQZ(pair_index)*m2y+BQY(pair_index))*m2z-BQX(pair_index)
         if (a_null(j,kkk).and.0.le.id.and.id.lt.max) then
          DO NUM1=1,NB
#ifndef _COMPLEX_
             SUMA(1:nb,num1)=SUMA(1:nb,num1)+a_(1:nb,j,kkk)*b_(num1,id)
#else
             SUMA(1:nb,num1)=SUMA(1:nb,num1)+conjg(a_(1:nb,j,kkk))*b_(num1,id)
#endif
           ENDDO
        endif
        ENDDO

 
        overlap(1:nb,1:nb,pair_index) = overlap(1:nb,1:nb,pair_index) + SUMA
        write(unit_mmn,235)KP(pair_index),KPB(pair_index),BQX1(pair_index),BQY1(pair_index), &
                   BQZ1(pair_index)

        DO NUM2=1,NB
         DO NUM1=1,NB
          write(unit_mmn,236) conjg(overlap(num2,num1,pair_index))
         ENDDO
        ENDDO
 4     CONTINUE
       call cputim(tt3)
       t_b=t_b+tt2-tt1
       t_prod=t_prod+tt3-tt2
       if (k1.eq.1) write(unit_out,*)'PLANEW TIME PER k-point:',t_b+t_prod
 3    CONTINUE 
      write(unit_out,*)'t_b=',t_b
      write(unit_out,*)'t_prod=',t_prod


 235  format(5i8)
 236  format(32f18.12)
      return

      !release the memory
      deallocate(FAC)
      deallocate(A_)
      deallocate(B_)
      deallocate(a_null)     
      deallocate(a_null_fin)

      END                                                               


!!/---
!! Local Variables:
!! mode: f90
!! End:
!!\---
!!
!! Time-stamp: <2015-06-30 21:37:47 assman@faepop23.tu-graz.ac.at>
