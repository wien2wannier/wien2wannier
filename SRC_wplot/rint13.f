!!! wien2wannier/SRC_wplot/rint13.f
!!!
!!! $Id: rint13.f 385 2015-06-01 13:08:18Z assmann $

      SUBROUTINE RINT13(REL,A,B,X,Y,S,JATOM)
      use radgrd                     
      use param
!     last changes: 01.11.00 ub (updating comments)
!
!     PERFORM RADIAL INTEGRALS
!     Int(0,Rmt) u(r) v(r) + 1/c^2 us(r) vs(r) r^2 dr
!
!     with c = 137.037 (274.074) in Hartree (Rydberg) units
!
!     For non-relativistic calculations c := 10^+11 Rydberg units.
!
!----------------------------------------------------------------------------
! Input:
! REL    .TRUE. for scalar relativistc calculations
! A(:)   r * u (r) on the radial mesh
! B(:)   r * us(r) on the radial mesh
! X(:)   r * v (r) on the radial mesh
! Y(:)   r * vs(r) on the radial mesh
! JATOM  the current type of atom
!
! from COMMON /RADGRD/
! RNOT(JATOM)  first radial mesh point
! DX  (JATOM)  logaritmic increment of the radial mesh
! JRI (JATOM)  number of radial mesh points
!
! Output:
! S    the value of the radial integral
!----------------------------------------------------------------------------
!
      IMPLICIT REAL(R8) (A-H,O-Z)
!
      LOGICAL REL
      DIMENSION A(NRAD),B(NRAD),X(NRAD),Y(NRAD)
!
!      COMMON /RADGRD/ RM(NRAD,NATO),RNOT(NATO),DX(NATO),JRI(NATO)
!
      CIN=1.331258D-5
      IF(.NOT.REL) CIN=1E-22
!     CONVERT FROM RYDBERG TO HARTREE UNITS
      CIN=CIN*4
!
      D=EXP(DX(JATOM))
      J=3-MOD(JRI(JATOM),2)
      J1=J-1
      R=RNOT(JATOM)*(D**(J-1))
      R1=R/D
      Z4=0
      Z2=0
   10 Z4=Z4+R*(A(J)*X(J)+CIN*B(J)*Y(J))
      R=R*D
      J=J+1
      IF(J.GE.JRI(JATOM)) GOTO 20
      Z2=Z2+R*(A(J)*X(J)+CIN*B(J)*Y(J))
      R=R*D
      J=J+1
      GOTO 10
   20 P1=RNOT(JATOM)*(A(1)*X(1)+CIN*B(1)*Y(1))
      P2=R1*(A(J1)*X(J1)+CIN*B(J1)*Y(J1))
      S=2*Z2+4*Z4+R*(A(J)*X(J)+CIN*B(J)*Y(J))+P2
      S=(DX(JATOM)*S+P1)/3.0D0
      IF(J1.GT.1) S=S+0.5D0*DX(JATOM)*(P1+P2)
      RETURN
      END


!!/---
!! Local Variables:
!! mode: f90
!! End:
!!\---
!!
!! Time-stamp: <2015-05-23 19:58:48 elias>
