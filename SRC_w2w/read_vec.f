!!! wien2wannier/SRC_w2w/read_vec.f
!!!
!!!    Read a ‘vector’ file
!!!
!!! Copyright 2010-2012 Jan Kuneš, Philipp Wissgott
!!!           2013-2015 Elias Assmann

SUBROUTINE read_vec(NEMIN,NEMAX,kkk,maxx,maxy,maxz,efermi)
  use param
  use const,  only: R8, Ryd_eV
  USE struct, only: nat
  USE xa3,    only: XK,YK,ZK, GX,GY,GZ, vecsz, A

  IMPLICIT none

  integer, intent(in)    :: NEMIN, NEMAX
  integer, intent(inout) :: maxx, maxy, maxz
  integer, intent(inout) :: kkk
  real(R8),intent(in)    :: efermi

  integer :: i, j, n, NB, NE, num

  real(r8) :: E(1000)

  CHARACTER(10)    BNAME

  DO I=1,NAT
     READ(unit_vector)
     READ(unit_vector)
  ENDDO

  maxx=0; maxy=0; maxz=0
  kpoint: do
     READ(unit_vector,END=998) XK(kkk+1),YK(kkk+1),ZK(kkk+1),BNAME,N,NE
     KKK=KKK+1
     vecsz(kkk)=N
     READ(unit_vector) (GX(I,kkk),GY(I,kkk),GZ(I,kkk), I=1,N)
     DO I=1,N
        IF (abs(GX(I,kkk)).gt.maxx) maxx=abs(GX(I,kkk))
        IF (abs(GY(I,kkk)).gt.maxy) maxy=abs(GY(I,kkk))
        IF (abs(GZ(I,kkk)).gt.maxz) maxz=abs(GZ(I,kkk))
     ENDDO

     DO J=1,NE
        READ(unit_vector)NUM,E(NUM)
        if (NUM.ge.NEMIN.and.NUM.le.nemax) then
           READ(unit_vector)(A(I,NUM-NEMIN+1,kkk),I=1,N)
        else
           READ(unit_vector)
        endif
     ENDDO

     DO NUM=NEMIN,NEMAX
        NB=NUM-NEMIN+1
        write(unit_eig,"(2I12,F22.16)")NB,kkk, ryd_ev*(E(NUM)-efermi)
     ENDDO
  end do kpoint

998 write(unit_out,*)'vector read in',kkk
END SUBROUTINE read_vec


!!/---
!! Local Variables:
!! mode: f90
!! End:
!!\---
!!
!! Time-stamp: <2016-07-05 14:20:37 assman@faepop71.tu-graz.ac.at>
