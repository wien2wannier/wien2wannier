!!! wien2wannier/SRC_w2w/read_vec.f
!!! 
!!!    Read a ‘vector’ file
!!!
!!! Copyright 2010-2012 Jan Kuneš, Philipp Wissgott
!!!                2013 Elias Assmann
!!!
!!! $Id: read_vec.f 191 2014-03-03 09:16:46Z assmann $

SUBROUTINE read_vec(NEMIN,NEMAX,kkk,maxx,maxy,maxz,efermi)
  use param
  use const,  only: R8, Ryd_eV
  USE struct, only: nat
  USE xa3,    only: XK,YK,ZK, GX,GY,GZ, size, A

  IMPLICIT none

  integer, intent(in)    :: NEMIN, NEMAX
  integer, intent(inout) :: maxx, maxy, maxz
  integer, intent(inout) :: kkk
  real(R8),intent(in)    :: efermi

  integer :: i, j, n, NB, NE, num

  real(r8) :: E(1000), Emist
  
  CHARACTER(10)    BNAME
  real(r8), allocatable ::  norm(:)
  real(r8), allocatable :: cmuss(:)

  DO I=1,NAT
     READ(unit_vector) EMIST
     READ(unit_vector) EMIST
  ENDDO

  kpoint: do
     READ(unit_vector,END=998) XK(kkk+1),YK(kkk+1),ZK(kkk+1),BNAME,N,NE
     KKK=KKK+1
     size(kkk)=N
     READ(unit_vector) (GX(I,kkk),GY(I,kkk),GZ(I,kkk), I=1,N)
     DO I=1,N
        IF (abs(GX(I,kkk)).gt.maxx) maxx=abs(GX(I,kkk))
        IF (abs(GY(I,kkk)).gt.maxy) maxy=abs(GY(I,kkk))
        IF (abs(GZ(I,kkk)).gt.maxz) maxz=abs(GZ(I,kkk))
     ENDDO

     allocate(cmuss(N))
     allocate(norm(NE))
     norm(:) = 0d0

     DO J=1,NE
        READ(unit_vector)NUM,E(NUM)
        if (NUM.ge.NEMIN.and.NUM.le.nemax) then
           READ(unit_vector)(A(I,NUM-NEMIN+1,kkk),I=1,N)
           !      do jb=1,N

           !       norm(J) = norm(J) + dabs(A(jb,NUM-NEMIN+1,kkk))**2
           !write(*,*)dabs(A(jb,NUM-NEMIN+1,kkk))**2,norm(J)
           !      enddo
        else
           READ(unit_vector)(cmuss(I),I=1,N)
           !     do jb=1,N
           !           norm(J) = norm(J) + dabs(cmuss(jb))**2
           !write(*,*)dabs(cmuss(jb))**2,norm(J)
           !     enddo

        endif
     ENDDO
     
     DO NUM=NEMIN,NEMAX
        NB=NUM-NEMIN+1
        write(unit_eig,"(2I12,F22.16)")NB,kkk, ryd_ev*(E(NUM)-efermi)
     ENDDO
     deallocate(cmuss)
     deallocate(norm)
  end do kpoint
  
998 write(unit_out,*)'vector read in',kkk
END SUBROUTINE read_vec

!!/---
!! Local Variables:
!! mode: f90
!! End:
!!\---
