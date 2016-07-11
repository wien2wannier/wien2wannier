!!! wien2wannier/SRC_w2w/outwin.f

subroutine OUTWIN(REL,V,RNOT,DH,JRI,EH,FL,VAL,SLO,Nodes,Z)
  !         Integration der skalarrel. Schroedingergleichung
  !
  !    Rydberg Einheiten

  use param, only: clight, Nrad
  use const, only: R8
  use uhelp, only: A, B

  implicit none

  real(R8) :: EH, fl, Z, V(Nrad), Rnot, dh, val, slo
  integer  :: jri, nodes
  logical  :: rel

  !  Input:
  !    EH    Energie in Hartree
  !    FL    Drehimpuls
  !    Z     Kernladung
  !    V     rad.sym. Potential in Hartree
  !    RNOT  erster radialer Netzpunkt
  !    DH    log. Schrittweite
  !    JRI   Anzahl radialer Netzpunkte
  intent(in)  :: EH, fl, Z, V, Rnot, dh, jri

  !  Output:
  !    VAL,SLO:  Wellenfunktion und Steigung am Kugelrand
  !    Nodes:    Anzahl Knoten
  intent(out) :: val, slo, nodes

  real(R8) :: D(2,3), Rnet(Nrad), C, E
  real(R8) :: zz, fllp1, s, sf, f0, aa, r, drdi, dg1,dg2,dg3, df1, df2, df3
  real(R8) :: phi, u, x, y, det, b1,b2
  integer  :: iiij, k

  real(R8), parameter :: H83 = 8/3._R8
  real(R8), parameter :: R83SQ = 64/9._R8
  real(R8), parameter :: R1 = 1/9._R8
  real(R8), parameter :: R2 = -5*R1
  real(R8), parameter :: R3 = 19*R1
  real(R8), parameter :: G0 = 1

  ! Hartree in Ryd
  E = 2*EH

  do iiij=1,JRI
     RNET(iiij)=RNOT*(exp(DH*(iiij-1)))
  end do

  Nodes = 0
  ZZ = Z + Z

  C = merge(2*clight, 1e10._R8, rel)

  FLLP1 = FL*(FL + 1)

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
  do K = 1,3
     R = RNET(K)
     DRDI = DH*R
     A(K) = (R**S)*G0
     B(K) = (R**SF)*F0
     D(1,K) = DRDI*A(K)*S/R
     D(2,K) = DRDI*B(K)*SF/R
  end do

  DG1 = D(1,1)
  DG2 = D(1,2)
  DG3 = D(1,3)
  DF1 = D(2,1)
  DF2 = D(2,2)
  DF3 = D(2,3)
  do K = 4, JRI
     R = RNET(K)
     DRDI = DH*R

     !       Faktor zwei vor V wegen Hartree-Rydberg !
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
  end do

  do iiij=1,JRI
     B(iiij)=B(iiij)*c/2.d0
  end do
  !
  VAL = A(JRI)/RNET(JRI)
  SLO = DG3/(DH*RNET(JRI))
  SLO = (SLO-VAL)/RNET(JRI)
end subroutine OUTWIN


!!/---
!! Local Variables:
!! mode: f90
!! End:
!!\---
!!
!! Time-stamp: <2016-07-11 12:03:07 assman@faepop71.tu-graz.ac.at>
