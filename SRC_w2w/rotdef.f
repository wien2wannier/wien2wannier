!!! wien2wannier/SRC_w2w/rotdef.f

SUBROUTINE ROTDEF
  !
  !     THE MATRICES  ROTIJ(3,3,INDEX)  TRANSFORM THE POSITION OF AN
  !     ATOM TO ITS CORRESPONDING POSITION OF AN EQUIVALENT ATOM.
  !
  use const,  only: R8
  use struct, only: Nat, mult, iord, iz, pos, trans, transij, rotij, lattic
  use clio,   only: croak
  use util,   only: string

  implicit none

  real(R8)    x(3),x1(3),toler,toler2,one
  INTEGER   i,m,index,index1,jatom,ncount
  DATA TOLER/1.D-7/,ONE/1.D0/
  toler2=1.5d0*toler
  INDEX=0
  NCOUNT=0
  DO JATOM=1,NAT
     INDEX1=INDEX+1
     DO M=1,MULT(JATOM)
        INDEX=INDEX+1
        DO I=1,IORD
           x(1:3)=0.d0
           x=MATMUL(TRANSPOSE(iz(1:3,1:3,i)),pos(1:3,index1))
           x(1:3)=x(1:3)+trans(1:3,i)
           x1(1:3)=MOD(ABS(X(1:3)-POS(1:3,INDEX))+toler,one)-toler
!           WRITE(*,*) 'JATOM,INDEX,I',JATOM,INDEX,I
!           WRITE(*,*) ABS(X1(1:3)),toler
           IF(MAXVAL(ABS(X1)) < TOLER2) THEN
              NCOUNT=NCOUNT+1
              TRANSIJ(1:3,INDEX)=TRANS(1:3,I)
              ROTIJ(1:3,1:3,INDEX)=IZ(1:3,1:3,I)
              GOTO 30
           END IF
           !....check positions for centered lattices
           if(lattic(1:1) == 'B') then
              x1(1:3)=mod(x1(1:3)+0.5d0+toler,one)-toler
              IF(MAXVAL(ABS(X1)) < TOLER2) THEN
                 NCOUNT=NCOUNT+1
                 TRANSIJ(1:3,INDEX)=TRANS(1:3,I)
                 ROTIJ(1:3,1:3,INDEX)=IZ(1:3,1:3,I)
                 GOTO 30
              END IF
           endif
           if(lattic(1:1) == 'F' .or. lattic(1:3) == 'CXY') then
              x1(1:2)=mod(x1(1:2)+0.5d0+toler,one)-toler
              IF(MAXVAL(ABS(X1)) < TOLER2) THEN
                 NCOUNT=NCOUNT+1
                 TRANSIJ(1:3,INDEX)=TRANS(1:3,I)
                 ROTIJ(1:3,1:3,INDEX)=IZ(1:3,1:3,I)
                 GOTO 30
              END IF
              x1(1:2)=mod(x1(1:2)+0.5d0,one)
           endif
           if(lattic(1:1) == 'F' .or. lattic(1:3) == 'CXZ') then
              x1(1)=mod(x1(1)+0.5d0+toler,one)-toler
              x1(3)=mod(x1(3)+0.5d0+toler,one)-toler
              IF(MAXVAL(ABS(X1)) < TOLER2) THEN
                 NCOUNT=NCOUNT+1
                 TRANSIJ(1:3,INDEX)=TRANS(1:3,I)
                 ROTIJ(1:3,1:3,INDEX)=IZ(1:3,1:3,I)
                 GOTO 30
              END IF
              x1(1)=mod(x1(1)+0.5d0,one)
              x1(3)=mod(x1(3)+0.5d0,one)
           endif
           if(lattic(1:1) == 'F' .or. lattic(1:3) == 'CYZ') then
              x1(2:3)=mod(x1(2:3)+0.5d0+toler,one)-toler
              IF(MAXVAL(ABS(X1)) < TOLER2) THEN
                 NCOUNT=NCOUNT+1
                 TRANSIJ(1:3,INDEX)=TRANS(1:3,I)
                 ROTIJ(1:3,1:3,INDEX)=IZ(1:3,1:3,I)
                 GOTO 30
              END IF
           end if
        ENDDO

        call croak("error in ROTDEF: no symmetry operation found &
             &for atoms "//trim(string(jatom))//" and "//trim(string(index))&
             &//", positions "//trim(string(POS(:, index1)))//" and " &
             &//trim(string(POS(:, INDEX))))
30      CONTINUE
     ENDDO
  ENDDO

  if (NCOUNT /= INDEX) call croak('error in ROTDEF: &
       &ROTIJ not defined for all atoms of basis.  NCOUNT = ' &
       &// trim(string(Ncount)))
END SUBROUTINE ROTDEF


!!/---
!! Local Variables:
!! mode: f90
!! End:
!!\---
!!
!! Time-stamp: <2016-07-15 16:34:06 assman@faepop71.tu-graz.ac.at>
