	SUBROUTINE read_vec(NEMIN,NEMAX,num_kpts,maxx,maxy,maxz)
  USE param
  USE struct
  USE xa3
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER *10    BNAME
      DIMENSION E(1000)
     real*8, allocatable ::  norm(:)
   integer jb
   real*8, allocatable :: cmuss(:)
  !_REAL      REAL*8 :: cdummy
  !_COMPLEX       COMPLEX*16 ::  cdummy

    MAXX=0
    MAXY=0
    MAXZ=0
    KKK=0
    
    DO I=1,NAT
     READ(9) EMIST
     READ(9) EMIST
    ENDDO
 4  CONTINUE
    KKK=KKK+1
    READ(9,END=998)XK(kkk),YK(kkk),ZK(kkk),BNAME,N,NE
    size(kkk)=N
    READ(9)(GX(I,kkk),GY(I,kkk),GZ(I,kkk),I=1,N)
    DO I=1,N
     IF (iabs(GX(I,kkk)).gt.maxx) maxx=iabs(GX(I,kkk))
     IF (iabs(GY(I,kkk)).gt.maxy) maxy=iabs(GY(I,kkk))
     IF (iabs(GZ(I,kkk)).gt.maxz) maxz=iabs(GZ(I,kkk))
    ENDDO
    
    allocate(cmuss(N))
    allocate(norm(NE))
    norm(:) = 0d0
    DO J=1,NE
     READ(9)NUM,E(NUM)
     if (NUM.ge.NEMIN.and.NUM.le.nemax) then
      READ(9)(A(I,NUM-NEMIN+1,kkk),I=1,N)
!      do jb=1,N

!       norm(J) = norm(J) + dabs(A(jb,NUM-NEMIN+1,kkk))**2
       !write(*,*)dabs(A(jb,NUM-NEMIN+1,kkk))**2,norm(J)
!      enddo
     else
     READ(9)(cmuss(I),I=1,N)
!     do jb=1,N
!           norm(J) = norm(J) + dabs(cmuss(jb))**2
           !write(*,*)dabs(cmuss(jb))**2,norm(J)
!     enddo
     
     endif
     
    ENDDO
    !write(*,*)'norm',norm
    DO NUM=NEMIN,NEMAX
     NB=NUM-NEMIN+1
     write(12,"(2I12,F22.16)")NB,kkk,13.6051413*E(NUM)
    ENDDO
    deallocate(cmuss)
    deallocate(norm)
    GOTO 4

 998  write(6,*)'vector read in',kkk-1,num_kpts
     !pw add 16.04.12
     maxx=maxx+1  
     maxy=maxy+1
     maxz=maxz+1
     !end pw add
     


    if(.not.(num_kpts.gt.(kkk-1))) stop 'inconsistent number of k-points'
    return
    END
    



