!!! wien2wannier/SRC_w2w/diracout.f

module     diracout_m; contains
subroutine diracout(rel,v,rnot,dstep,nmax,eh,nqk,val,slo,nodes,z)
!!! Integration of Dirac equation

  use param, only: clight, unit_out, Nrad
  use const, only: R8
  ! dp    =  large component of the solution of the dirac equation
  ! dq    =  small component of the solution
  use uhelp, only: dp => A, dq => B
  use PS1,   only: dep, deq, db, dvc, dsal, dk, dm

  !! procedure includes
  use inth_m
  use inouh_m

  implicit none

  !  Input:
  !    rel    switch for relativ. - nonrelativ. calculation
  !    V      rad.sym. potential in Hartree, V = potential*r
  !    Rnot   first radial meshpoint
  !    dstep  log. step
  !    Nmax   number of radial meshpoints
  !    EH     energy in hartree
  !    Nqk    relativistic quantum number kappa
  !    Z      charge of nucleus
  real(R8), intent(in)  :: V(Nrad), Rnot, dstep, EH, Z
  logical,  intent(in)  :: rel
  integer,  intent(in)  :: Nmax, Nqk

  !  Output:
  !    val,slo:  Wellenfunktion und Steigung am Kugelrand
  !    nodes:    nomber of nodes
  real(R8), intent(out) :: val
  integer,  intent(out) :: nodes

  real(R8), parameter :: dkoef=1 / 720._R8, test=1.e-8_R8
  real(R8), parameter :: atom_mass(103) = (/ &
       &   1.0_R8,   4.0_R8,   6.9_R8,   9.0_R8,  10.8_R8,  12.0_R8, &
       &  14.0_R8,  16.0_R8,  19.0_R8,  20.2_R8,  23.0_R8,  24.3_R8, &
       &  27.0_R8,  28.1_R8,  31.0_R8,  32.0_R8,  35.4_R8,  40.0_R8, &
       &  39.1_R8,  40.0_R8,  45.0_R8,  47.9_R8,  50.9_R8,  52.0_R8, &
       &  54.9_R8,  55.8_R8,  58.9_R8,  58.7_R8,  63.5_R8,  65.4_R8, &
       &  69.7_R8,  72.6_R8,  74.9_R8,  79.0_R8,  79.9_R8,  83.8_R8, &
       &  85.5_R8,  87.6_R8,  88.9_R8,  91.2_R8,  92.9_R8,  95.9_R8, &
       &  98.0_R8, 101.1_R8, 102.9_R8, 106.4_R8, 107.9_R8, 112.4_R8, &
       & 114.8_R8, 118.7_R8, 121.8_R8, 127.6_R8, 126.9_R8, 131.3_R8, &
       & 132.9_R8, 137.3_R8, 138.9_R8, 140.1_R8, 140.9_R8, 144.2_R8, &
       & 145.0_R8, 150.4_R8, 152.0_R8, 157.3_R8, 158.9_R8, 162.5_R8, &
       & 164.9_R8, 167.3_R8, 168.9_R8, 173.0_R8, 175.0_R8, 178.5_R8, &
       & 180.9_R8, 183.8_R8, 186.2_R8, 190.2_R8, 192.2_R8, 195.1_R8, &
       & 197.0_R8, 200.6_R8, 204.4_R8, 207.2_R8, 209.0_R8, 209.0_R8, &
       & 210.0_R8, 222.0_R8, 223.0_R8, 226.0_R8, 227.0_R8, 232.0_R8, &
       & 231.0_R8, 238.0_R8, 237.0_R8, 244.0_R8, 243.0_R8, 247.0_R8, &
       & 247.0_R8, 251.0_R8, 252.0_R8, 257.0_R8, 258.0_R8, 259.0_R8, &
       & 262.0_R8 /)

  !      DR   =    radial mesh
  real(R8) :: dv(nrad),  dr(NRAD)
  real(R8) :: d1, dfl, dq1, Rnuc, dval, slo
  integer  :: i, nuc

  !rschmid
  !   Set up radial mesh.
  !rschmid
  do i=1,nmax
     DR(i)=RNOT*(exp(DSTEP*(i-1.d0)))
  enddo

  if (rel) then
     dvc = clight
  else
     dvc = 1e10_R8
  endif
  dsal = 2*dvc
  db = eh/dvc
  dk = nqk
  dm=dstep*dkoef

  do i=1,nmax
     dv(i) = v(i)/dr(i)
  enddo
  !rschmid
  !  Behavior of the solution at the origin
  !rschmid

  !jk   finite size of the nucleus
  rnuc=2.2677e-05_R8*(atom_mass(int(z))**(1/3._R8))
  write(unit_out,*)'amass, r0:',atom_mass(int(z)),rnuc
  do i=1,nmax
     d1=rnot*exp(DSTEP*(i-1.d0))
     if (d1 >= rnuc) exit
  end do
  nuc=I
  write(unit_out,*)'nuc=',nuc
  if (nuc <= 0) then
     dfl = sqrt(nqk*nqk-z*z/(dvc*dvc))
  else
     dfl=nqk*nqk
     do i=1,nuc
        dv(i)=dv(i)+z/dr(i)+z*((dr(i)/dr(nuc))**2-3)/(2*dr(nuc))
     end do
  end if
  dq1 = nqk/iabs(nqk)

  !rschmid
  !  Determine expansion of the potential at the origin.
  !rschmid
  call inouh(dp,dq,dr,dq1,dfl,dv(1),Z,TEST,nuc)

  !rschmid
  !  Set up boundary conditions using a small r expansion.
  !rschmid

  nodes = 0
  do i=1,5
     dval=dr(i)**dfl
     if (i /= 1) then
        if (dp(i-1) /= 0) then
           if ((dp(i)/dp(i-1)) <= 0) then
              nodes=nodes+1
           endif
        endif
     endif
     dp(i) = dp(i)*dval
     dq(i) = dq(i)*dval
     dep(i)=dep(i)*dval
     deq(i)=deq(i)*dval
  enddo

  !rschmid
  !    Perform outward integration of dirac equation
  !rschmid

  do i = unit_out, nmax
     dp(i) = dp(i-1)
     dq(i) = dq(i-1)
     call inth (dp(i),dq(i),dv(i),dr(i))
     if (dp(i-1) /= 0) then
        if ((dp(i)/dp(i-1)) > 0) then
           nodes=nodes+1
        endif
     endif
  enddo


  val = dp(nmax)/dr(nmax)
  slo = dep(5)/(dstep*dr(nmax))/dvc*2.d0
  slo = (slo-val)/dr(nmax)
end subroutine diracout
end module     diracout_m

!!/---
!! Local Variables:
!! mode: f90
!! End:
!!\---
!!
!! Time-stamp: <2016-07-18 17:29:39 assman@faepop71.tu-graz.ac.at>
