!!! wien2wannier/SRC_wplot/ylm.f

subroutine YLM(V,LMAX,Y)
  use const, only: R8, C16, TAU

  implicit none

  integer,      intent(in)  :: LMAX
  real(R8),     intent(in)  :: V(3)
  complex(C16), intent(out) :: Y(*)

!     ..................................................................
! 1.     PROGRAM UNIT 'YLM'
!           Calculates spherical harmonics
!           FORTRAN 77 SUBROUTINE
!
! 2.     PURPOSE
!           The spherical harmonics (Condon and Shortley convention)
!             Y(0,0),Y(1,-1),Y(1,0),Y(1,1),Y(2,-2) ... Y(LMAX,LMAX)
!           for vector V (given in Cartesian coordinates)
!           are calculated. In the Condon Shortley convention the
!           spherical harmonics are defined as
!                             +------+
!                        m    |   1     m              im(Phi)
!           Y(l,m) = (-1)  -+ | -----  P (cos(Theta)) e
!                            \| 2(Pi)   l
!                  m
!           where P (cos(Theta)) is the normalized Associated Legendre
!                  l
!           function. Thus,
!                                          m      *
!                            Y(l,-m) = (-1) Y(l,m)
!
!
! 3.     USAGE
!           DOUBLE PRECISION V(3), Y(5*5)
!           V(1) = ...
!           V(2) = ...
!           V(3) = ...
!           CALL YLM(V,4,Y)
!
!        ARGUMENT-DESCRIPTION
!           V      - DOUBLE PRECISION vector, dimension 3        (input)
!                    Must be given in Cartesian coordinates.
!                    Conversion of V to polar coordinates gives the
!                    angles Theta and Phi necessary for the calculation
!                    of the spherical harmonics.
!           LMAX   - INTEGER value                               (input)
!                    upper bound of L for which spherical harmonics
!                    will be calculated
!                    constraint:
!                       LMAX .GE. 0 (not checked)
!           Y      - COMPLEX*16 array, dimension (LMAX+1)**2    (output)
!                    contains the calculated spherical harmonics
!                    Y(1)                   for L .EQ. 0 (M = 0)
!                    Y(2), ..., Y(4)        for L .EQ. 1 (M = -1, 0, 1)
!                    ...
!                    Y(LMAX*LMAX+1), ..., Y((LMAX+1)*(LMAX+1))
!                                           for L .EQ. LMAX
!                                              (M = -L,...,L)
!                    constraint:
!                       Dimension of Y .GE. (LMAX+1)**2 (not checked)
!
!        USED SUBROUTINES (DIRECTLY CALLED)
!           none
!
!        INDIRECTLY CALLED SUBROUTINES
!           none
!
!        UTILITY-SUBROUTINES (USE BEFOREHAND OR AFTERWARDS)
!           none
!
!        INPUT/OUTPUT (READ/WRITE)
!           none
!
!        MACHINENDEPENDENT PROGRAMPARTS
!           Type COMPLEX*16 is used which does not conform to the
!           FORTRAN 77 standard.
!           Also the non-standard type conversion function DCMPLX()
!           is used which combines two double precision values into
!           one double complex value.
!
! 4.     REMARKS
!           none
!
! 5.     METHOD
!           The basic algorithm used to calculate the spherical
!           harmonics for vector V is as follows:
!
!           Y(0,0)
!           Y(1,0)
!           Y(1,1)
!           Y(1,-1) = -Y(1,1)
!           DO L = 2, LMAX
!              Y(L,L)   = f(Y(L-1,L-1)) ... Formula 1
!              Y(L,L-1) = f(Y(L-1,L-1)) ... Formula 2
!              DO M = L-2, 0, -1
!                 Y(L,M) = f(Y(L-1,M),Y(L-2,M)) ... Formula 2
!                 Y(L,-M)= (-1)**M*Y(L,M)
!              ENDDO
!           ENDDO
!
!           In the following the necessary recursion formulas and
!           starting values are given:
!
!        Start:
!                        +------+
!                        |   1
!           Y(0,0) =  -+ | -----
!                       \| 4(Pi)
!
!                                   +------+
!                                   |   3
!           Y(1,0) =  cos(Theta) -+ | -----
!                                  \| 4(Pi)
!
!                                     +------+
!                                     |   3    i(Phi)
!           Y(1,1) =  - sin(Theta) -+ | ----- e
!                                    \| 8(Pi)
!
!        Formula 1:
!
!           Y(l,l) =
!                           +--------+
!                           | (2l+1)   i(Phi)
!            -sin(Theta) -+ | ------  e       Y(l-1,l-1)
!                          \|   2l
!
!        Formula 2:
!                                  +---------------+
!                                  |  (2l-1)(2l+1)
!           Y(l,m) = cos(Theta) -+ | -------------- Y(l-1,m)  -
!                                 \|   (l-m)(l+m)
!
!                                    +--------------------+
!                                    |(l-1+m)(l-1-m)(2l+1)
!                              -  -+ |-------------------- Y(l-2,m)
!                                   \|  (2l-3)(l-m)(l+m)
!
!        Formula 3: (not used in the algorithm because of the division
!                    by sin(Theta) which may be zero)
!
!                                    +--------------+
!                      cos(Theta)    |  4(m+1)(m+1)   -i(Phi)
!           Y(l,m) = - ---------- -+ | ------------  e       Y(l,m+1) -
!                      sin(Theta)   \| (l+m+1)(l-m)
!
!                                    +--------------+
!                                    |(l-m-1)(l+m+2)  -2i(Phi)
!                              -  -+ |-------------- e        Y(l,m+2)
!                                   \| (l-m)(l+m+1)
!
! 6.     DATE
!           26. April 1994                                   Version 1.2
!
!        INSTITUT FUER TECHNISCHE ELEKTROCHEMIE            --  TU VIENNA
!     ..................................................................

  integer  :: I2L, I4L2, INDEX, INDEX2, L, M, MSIGN
  real(R8) :: A, B, C, AB, ABC, ABMAX, ABCMAX, D4LL1C, D2L13
  real(R8) :: COSTH, SINTH, COSPH, SINPH, TEMP1, TEMP2, TEMP3
  real(R8) :: YLLR, YLL1R, YL1L1R, YLMR, YLLI, YLL1I, YL1L1I, YLMI

! Y(0,0)
  YLLR = 1/sqrt(2*TAU)
  YLLI = 0
  Y(1) = cmplx(YLLR, YLLI, R8)

! continue only if spherical harmonics for (L .GT. 0) are desired
  if (LMAX <= 0) return

! calculate sin(Phi), cos(Phi), sin(Theta), cos(Theta)
! Theta, Phi ... polar angles of vector V
  ABMAX  = max(abs(V(1)),abs(V(2)))
  if (ABMAX > 0) then
     A = V(1)/ABMAX
     B = V(2)/ABMAX
     AB = sqrt(A*A+B*B)
     COSPH = A/AB
     SINPH = B/AB
  else
     COSPH = 1
     SINPH = 0
  endif
  ABCMAX = max(ABMAX,abs(V(3)))
  if (ABCMAX > 0) then
     A = V(1)/ABCMAX
     B = V(2)/ABCMAX
     C = V(3)/ABCMAX
     AB = A*A + B*B
     ABC = sqrt(AB + C*C)
     COSTH = C/ABC
     SINTH = sqrt(AB)/ABC
  else
     COSTH = 1
     SINTH = 0
  endif

! Y(1,0)
  Y(3) = cmplx(sqrt(3._R8)*YLLR*COSTH, 0, R8)

! Y(1,1) ( = -DCONJG(Y(1,-1)))
  TEMP1 = -sqrt(1.5D+0)*YLLR*SINTH
  Y(4) = DCMPLX(TEMP1*COSPH,TEMP1*SINPH)
  Y(2) = -DCONJG(Y(4))

  do L = 2, LMAX
     INDEX  = L*L+1
     INDEX2 = INDEX + 2*L
     MSIGN  = 1 - 2*mod(L,2)

!    YLL = Y(L,L) = f(Y(L-1,L-1)) ... Formula 1
     YL1L1R = dble(Y(INDEX-1))
     YL1L1I = aimag(Y(INDEX-1))
     TEMP1 = -sqrt(dble(2*L+1)/dble(2*L))*SINTH
     YLLR = TEMP1*(COSPH*YL1L1R - SINPH*YL1L1I)
     YLLI = TEMP1*(COSPH*YL1L1I + SINPH*YL1L1R)
     Y(INDEX2) = DCMPLX(YLLR,YLLI)
     Y(INDEX)  = MSIGN*DCONJG(Y(INDEX2))
     INDEX2 = INDEX2 - 1
     INDEX  = INDEX  + 1

!    YLL1 = Y(L,L-1) = f(Y(L-1,L-1)) ... Formula 2
!           (the coefficient for Y(L-2,L-1) in Formula 2 is zero)
     TEMP2 = sqrt(dble(2*L+1))*COSTH
     YLL1R = TEMP2*YL1L1R
     YLL1I = TEMP2*YL1L1I
     Y(INDEX2) = DCMPLX(YLL1R,YLL1I)
     Y(INDEX)  = -MSIGN*DCONJG(Y(INDEX2))
     INDEX2 = INDEX2 - 1
     INDEX  = INDEX  + 1

     I4L2 = INDEX2 - 4*L + 2
     I2L  = INDEX2 - 2*L
     D4LL1C = COSTH*sqrt(dble(4*L*L-1))
     D2L13  = -sqrt(dble(2*L+1)/dble(2*L-3))

     do M = L - 2, 0, -1
!       YLM = Y(L,M) = f(Y(L-2,M),Y(L-1,M)) ... Formula 2
        TEMP1 = 1.0D+0/sqrt(dble((L+M)*(L-M)))
        TEMP2 = D4LL1C*TEMP1
        TEMP3 = D2L13*sqrt(dble((L+M-1)*(L-M-1)))*TEMP1
        YLMR = TEMP2*dble(Y(I2L))  + TEMP3*dble(Y(I4L2))
        YLMI = TEMP2*aimag(Y(I2L)) + TEMP3*aimag(Y(I4L2))
        Y(INDEX2) = DCMPLX(YLMR,YLMI)
        Y(INDEX)  = MSIGN*DCONJG(Y(INDEX2))

        MSIGN  = -MSIGN
        INDEX2 = INDEX2 - 1
        INDEX  = INDEX  + 1
        I4L2   = I4L2   - 1
        I2L    = I2L    - 1
     end do
  end do
end subroutine YLM


!!/---
!! Local Variables:
!! mode: f90
!! End:
!!\---
!!
!! Time-stamp: <2016-07-15 14:14:24 assman@faepop71.tu-graz.ac.at>
