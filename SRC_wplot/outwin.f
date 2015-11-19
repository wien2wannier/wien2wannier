!!! wien2wannier/SRC_wplot/outwin.f
!!!
!!! $Id: outwin.f 166 2014-02-03 09:39:24Z assmann $

SUBROUTINE OUTWIN(REL,V,RI,DH,JRI,EH,FL,VAL,SLO,Nodes,Z) 
!     Integration of the scalar-relativistic Schroedinger equation
!     with psi(r) = u_l(|r|) Y_lm(r/|r|)
! ----------------------------------------------------------------
!  Input:
!    REL   .TRUE. for skalarrelativistic calculation
!    V(:)  radialsymmetric potential in Hartree
!    RI(:) radial mesh points
!    DH    (logical) step width
!    JRI   number of radial mesh points
!    EH    energy in Hartree 
!    FL    angular momentum
!    Z     charge of nucleus
!
!  Output:
!    VAL       wave function u_l(r) at MT sphere r = Rmt
!    SLO       descent d/dr u_l(r) at MT sphere r = Rmt
!    Nodes     number of nodes
!
!  COMMON /WORK/
!    A(:)   r * u _l(r)     at mesh points
!    B(:)   r * us_l(r) * c at mesh points (2nd rel. component)
!
!           Note, that here <(u,us)|(v,vs)> := <u|v> + <us|vs>  
!
!  Rydberg units
!
! ----------------------------------------------------------------
  use const

  use param
  
  IMPLICIT REAL(R8) (A-H,O-Z)

  LOGICAL       REL
  DIMENSION     V(NRAD),RI(NRAD),D(2,3)
  COMMON /WORK1/ A(NRAD),B(NRAD)

! Hartree in Ryd
  E=EH*2.d0

  Nodes = 0
  ZZ = Z + Z
  C = 274.074D0
  if(.not.rel) C=1.D+10

  FLLP1 = FL*(FL + 1.d0)
  R83SQ = 64.D0/9.D0
  R1 = 1.D0/9.D0
  R2 = -5.D0*R1
  R3 = 19.D0*R1
  H83 = 8.D0/3.D0

  G0 = 1.d0
  IF (Z .LT. 0.9D0) THEN
     S = FL+1.d0
     SF = FL
     F0 = FL/C
  ELSE
     AA = ZZ/C
     S = DSQRT(FLLP1 + 1.D0 - AA*AA)
     SF = S
     F0 = G0*(S - 1.D0)/AA
  ENDIF
  DO K = 1,3
     R = RI(K)
     DRDI = DH*R
     A(K) = (R**S)*G0
     B(K) = (R**SF)*F0
     D(1,K) = DRDI*A(K)*S/R
     D(2,K) = DRDI*B(K)*SF/R
  END DO

  DG1 = D(1,1)
  DG2 = D(1,2)
  DG3 = D(1,3)
  DF1 = D(2,1)
  DF2 = D(2,2)
  DF3 = D(2,3)
  DO K = 4, JRI
     R = RI(K)
     DRDI = DH*R

!    factor of 2 before V because of Hartree-Rydberg !
     PHI = (E - 2.d0*V(K)/R)*DRDI/C
     U = DRDI*C + PHI
     X = -DRDI/R
     Y = -FLLP1*X*X/U + PHI
     DET = R83SQ - X*X + U*Y
     B1 = A(K-1)*H83 + R1*DG1 + R2*DG2 + R3*DG3
     B2 = B(K-1)*H83 + R1*DF1 + R2*DF2 + R3*DF3
     A(K) = (B1*(H83-X) + B2*U)/DET
     B(K) = (B2*(H83+X) - B1*Y)/DET
     IF (A(K)*A(K-1) .LT. 0D0) Nodes = Nodes + 1
     DG1 = DG2
     DG2 = DG3
     DG3 = U*B(K) - X*A(K)
     DF1 = DF2
     DF2 = DF3
     DF3 = X*B(K) - Y*A(K)
  END DO

  DO iiij=1,JRI
     B(iiij)=B(iiij)*c/2.d0
  END DO

  VAL = A(JRI)/RI(JRI)
  SLO = DG3/(DH*RI(JRI))
  SLO = (SLO-VAL)/RI(JRI) 
END SUBROUTINE OUTWIN

!!/---
!! Local Variables:
!! mode: f90
!! End:
!!\---
