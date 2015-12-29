!!! wien2wannier/SRC_wplot/rint13.f

      subroutine RINT13(rel, A, B, X, Y, S, jatom, stru)
      use radgrd,    only: dx
      use param,     only: DPk, Nrad, clight
      use structmod, only: struct_t

      implicit none

      logical,        intent(in)  :: rel
      real(DPk),      intent(in)  :: A(Nrad), B(Nrad), X(Nrad), Y(Nrad)
      real(DPk),      intent(out) :: S
      integer,        intent(in)  :: jatom
      type(struct_t), intent(in)  :: stru

      real(DPk) :: cin, d, R, R1, Z2, Z4, P1, P2
      integer   :: j, j1

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
! from RADGRD
! DX  (JATOM)  logaritmic increment of the radial mesh
!
! Output:
! S    the value of the radial integral
!----------------------------------------------------------------------------

      CIN = 1/clight**2
      IF(.NOT.REL) CIN=4E-22    ! legacy value
!     CONVERT FROM RYDBERG TO HARTREE UNITS
!
      D=EXP(DX(JATOM))
      J=3-MOD(STRU%NPT(JATOM),2)
      J1=J-1
      R=STRU%R0(JATOM)*(D**(J-1))
      R1=R/D
      Z4=0
      Z2=0
   10 Z4=Z4+R*(A(J)*X(J)+CIN*B(J)*Y(J))
      R=R*D
      J=J+1
      IF(J.GE.STRU%NPT(JATOM)) GOTO 20
      Z2=Z2+R*(A(J)*X(J)+CIN*B(J)*Y(J))
      R=R*D
      J=J+1
      GOTO 10
   20 P1=STRU%R0(JATOM)*(A(1)*X(1)+CIN*B(1)*Y(1))
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
!! Time-stamp: <2015-12-29 19:46:44 assman@faepop36.tu-graz.ac.at>
