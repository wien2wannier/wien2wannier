!!! wien2wannier/SRC_w2w/diracout.f

subroutine diracout(rel,v,rnot,dstep,nmax,eh,nqk,val,slo,nodes,z)

! Integration of Dirac equation
! -----------------------------
      use param, only: clight, unit_out, Nrad
      use const, only: R8
!     dp    =  large component of the solution of the dirac equation
!     dq    =  small component of the solution
      use uhelp, only: dp => A, dq => B
      use PS1,   only: dep, deq, db, dvc, dsal, dk, dm

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
           &   1.0,   4.0,   6.9,   9.0,  10.8,  12.0,  14.0,  16.0,  19.0, &
           &  20.2,  23.0,  24.3,  27.0,  28.1,  31.0,  32.0,  35.4,  40.0, &
           &  39.1,  40.0,  45.0,  47.9,  50.9,  52.0,  54.9,  55.8,  58.9, &
           &  58.7,  63.5,  65.4,  69.7,  72.6,  74.9,  79.0,  79.9,  83.8, &
           &  85.5,  87.6,  88.9,  91.2,  92.9,  95.9,  98.0, 101.1, 102.9, &
           & 106.4, 107.9, 112.4, 114.8, 118.7, 121.8, 127.6, 126.9, 131.3, &
           & 132.9, 137.3, 138.9, 140.1, 140.9, 144.2, 145.0, 150.4, 152.0, &
           & 157.3, 158.9, 162.5, 164.9, 167.3, 168.9, 173.0, 175.0, 178.5, &
           & 180.9, 183.8, 186.2, 190.2, 192.2, 195.1, 197.0, 200.6, 204.4, &
           & 207.2, 209.0, 209.0, 210.0, 222.0, 223.0, 226.0, 227.0, 232.0, &
           & 231.0, 238.0, 237.0, 244.0, 243.0, 247.0, 247.0, 251.0, 252.0, &
           & 257.0, 258.0, 259.0, 262.0 /)

!      DR   =    radial mesh
      real(R8) :: dv(nrad),  dr(NRAD)
      real(R8) :: d1, dfl, dq1, Rnuc, dval, slo
      integer  :: i, nuc, nstop

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
      do 30 i=1,nuc
         dv(i)=dv(i)+z/dr(i)+z*((dr(i)/dr(nuc))**2-3.)/(2*dr(nuc))
30    continue
      end if
      dq1 = nqk/iabs(nqk)

!rschmid
!  Determine expansion of the potential at the origin.
!rschmid
      CALL INOUH (dp,dq,dr,dq1,dfl,dv(1),Z,TEST,nuc,NSTOP)

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

      do i = unit_out, nmax
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
!! Time-stamp: <2016-07-06 13:02:26 assman@faepop71.tu-graz.ac.at>
