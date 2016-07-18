!!! wien2wannier/SRC_w2w/gaunt1.f

module                 gaunt1_m; contains
real(R8) pure function gaunt1(LP,L,LS,MP,M,MS)
!
  use assleg, only: YR, N
  use const,  only: R8

  implicit none
!
!        Arguments
!
  integer, intent(in) :: L, LP, LS, M, MP, MS
!
!     ..................................................................
!
!        GAUNT computes the integral of
!           CONJG(Y(LP,MP))*Y(L,M)*Y(LS,MS) for LP+L+LS .LE. 23
!        using gaussian quadrature with N=12 as given by
!           M. Abramowitz and I.A. Stegun,
!           'Handbook of Mathematical Functions',
!           NBS Applied Mathematics Series 55 (1968), pages 887 and 916
!        written by Bruce Harmon based on suggestion by W. Rudge
!        Iowa State Sept.1973
!
!        extended by M. Weinert and E. Wimmer
!        Northwestern University March 1980
!
!     ..................................................................
!
!      INTEGER            MAXDIM, N
!      PARAMETER          (MAXDIM = 81, N = 6)
!
!        Local Scalars
!
      integer  :: I, IL, ILP, ILS
      real(R8) :: S
      real(R8), parameter :: W(N) = (/ &
           0.24914704581340D+0, 0.23349253653836D+0, &
           0.20316742672307D+0, 0.16007832854335D+0, &
           0.10693932599532D+0, 0.04717533638651D+0 /)

      IL = L*(L+1) + M + 1
      ILP = LP*(LP+1) + MP + 1
      ILS = LS*(LS+1) + MS + 1
      S = 0
      do I = 1, N
         S = S + W(I)*YR(I,ILP)*YR(I,IL)*YR(I,ILS)
      end do
      GAUNT1 = S
!
end function gaunt1
end module   gaunt1_m


!!/---
!! Local Variables:
!! mode: f90
!! End:
!!\---
!!
!! Time-stamp: <2016-07-18 13:56:26 assman@faepop71.tu-graz.ac.at>
