!!! wien2wannier/SRC_wplot/diracout.f

subroutine diracout(rel,v,rnot,dstep,nmax,eh,nqk,val,slo,nodes,z)
!rschmid
!         Integration of Dirac equation.
! 
!  Input:
! 
!    rel    switch for relativ. - nonrelativ. calculation
!    v      rad.sym. potential in Hartree
!    rnot   first radial meshpoint
!    dstep  log. step
!    nmax    number of radial meshpoints
!    eh     energy in hartree
!    nqk    relativistic quantum number kappa 
!    z      charge of nucleus
!
!  Output: 
!
!    val,slo:  Wellenfunktion und Steigung am Kugelrand
!    nodes:    nomber of nodes
!
!rschmid
!
! ----------------------------------------------------------------
  use ams,   only: atom_mass
  use ps1,   only: dep, deq, db=>dd, dvc, dsal, dk, dm=>dm1
  use param, only: DPk, Nrad, unit_out, clight

  implicit none

  logical :: rel

  real(DPk) :: v(Nrad), rnot, dstep, eh, val, slo, z
  integer   :: nmax, nqk, nodes

  intent(in)  :: rel, v, rnot, dstep, nmax, eh, nqk, z
  intent(out) :: val, slo, nodes

!rschmid
!     V     =   potential*r
!      DR   =    radial mesh
!     dp    =  large component of the solution of the dirac equation
!     dq    =  small component of the solution 
!rschmid
  real(DPk), dimension(Nrad) :: dr, dp, dq, dv

  save :: dp, dq

  real(DPk), parameter :: DKOEF = 1/720._DPk

! DEP,DEQ DERIVEES DE DP ET DQ   DB=ENERGIE/DVC    DVC VITESSE DE LA
! LUMIERE EN U.A.   DSAL=2.*DVC   DK NOMBRE QUANTIQUE KAPPA
! DM=PAS EXPONENTIEL/720., DKOEF=1./720.

!rschmid
!  The name of dm should be changed to avoid a name collision
!  with dm in inouh
!rschmid

  integer   :: i, nuc
  real(DPk) :: rnuc, d1, test, dval, dfl, dq1

!rschmid
!   Set up radial mesh.
!rschmid
      do i=1,nmax
         DR(i)=RNOT*(exp(DSTEP*(i-1.d0)))
      enddo

      if (rel) then
         dvc = clight
      else 
         dvc = 1.e10_DPk
      endif

      dsal = 2.d0*dvc
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
      rnuc=2.2677e-05_DPk * (atom_mass(int(z))**(1/3.))
      write(unit_out,*)'amass, r0:',atom_mass(int(z)),rnuc
      do 10 i=1,nmax
       d1=rnot*exp(DSTEP*(i-1.d0))
      if (d1.ge.rnuc) goto 20
 10   continue
 20   nuc=I
      write(unit_out,*)'nuc=',nuc
      if (nuc.le.0) then
      dfl = sqrt(nqk*nqk-z*z/(dvc*dvc))
      else
      dfl=nqk*nqk
      do i=1,nuc
         dv(i)=dv(i)+z/dr(i)+z*((dr(i)/dr(nuc))**2-3.)/(2*dr(nuc))
      end do
      end if
      dq1 = nqk/iabs(nqk)
     
!rschmid
!  Determine expansion of the potential at the origin.
!rschmid
      test =1.e-8
    
      CALL INOUH (dp,dq,dr,dq1,dfl,dv(1),Z,TEST,nuc)

!rschmid
!  Set up boundary conditions using a small r expansion.
!rschmid

      nodes = 0
      do i=1,5
        dval=dr(i)**dfl
        if (i.ne.1) then
          if (dp(i-1).ne.0.) then
             if ((dp(i)/dp(i-1)).le.0.) then
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

      do i = 6, nmax
        dp(i) = dp(i-1)
        dq(i) = dq(i-1)
        call inth (dp(i),dq(i),dv(i),dr(i))
        if (dp(i-1).ne.0.) then
          if ((dp(i)/dp(i-1)).gt.0.) then
            nodes=nodes+1
          endif
        endif
      enddo


      val = dp(nmax)/dr(nmax)
      slo = dep(5)/(dstep*dr(nmax))/dvc*2.d0
      slo = (slo-val)/dr(nmax) 

      RETURN
      END


!!/---
!! Local Variables:
!! mode: f90
!! End:
!!\---
!!
!! Time-stamp: <2015-12-29 19:15:18 assman@faepop36.tu-graz.ac.at>
