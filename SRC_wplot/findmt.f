!!! wien2wannier/SRC_wplot/findmt.f

subroutine FINDMT(P, atms, stru, iAt, iLat, iR, R)
      use struct,    only: RMT, POS
      use radgrd,    only: dx
      use const,     only: DPk
      use latt,      only: br2
      use structmod, only: struct_t

      implicit none

      type(struct_t), intent(in)  :: stru
      real(DPk),      intent(in)  :: P(3), atms(3, stru%Nneq)
      integer,        intent(out) :: iAt, iLat(3), iR
      real(DPk),      intent(out) :: R(3)

      integer   :: i, j, jatom, jx, jy, jz, ilow, iup
      real(DPk) :: RR, R2, RMT2, T

!
! checks whether a given real space point r falls into a muffin tin
! sphere around R + R_a with a running over all atoms in the unit cell
! --------------------------------------------------------------------
! Input:
! P(:)  -- the real space point in primitive fractional coordinates
! NAT   -- total number of atoms per unit cell
! ATMS  -- the size of the smallest primitive unit cells surrounding
!          the muffin tin spheres of each group of sym.-eq. atoms
!
! module LATT    data on the Bravais lattice
!
! Output:
! IAT   -- the atom a the muffin tin sphere belongs to
!          or zero if the given point is in the interstitial region
! ILAT  -- the lattice displacement R = Sum(i) N_i a_i in primitive
!          fractional coordinates [ undefined if IAT = 0 ]
! IR    -- radial mesh index i such that |r - R - R_a| in [r_i,r_i+1]
!          [ undefined if IAT = 0 ]
! R(:)  -- if in interstitial:
!          the absolute position r in Cartesian coordinates
!          if in muffin tin spheres:
!          the relative positition r - R - R_a in Cartesian corrdinates
!
!          Here a_i are the primitive real space lattice vectors !
! --------------------------------------------------------------------

      DIMENSION T(3),ILOW(3),IUP(3)
!
      DO 10 IAT=1,STRU%NAT
        JATOM = stru%neq2at(iat)
        RMT2 = RMT(JATOM)*RMT(JATOM)
!
!       << setup search area for atomic sphere at R0+R >>
! --------------------------------------------------------------------
! Letting R = Sum(i) N_i a_i the search area for R is given by
! all integers N_i in [ x_i - R0_i - s_i , x_i - R0_i + s_i ]
!
! to account for numerical noise let s_i = s_i + 0.0001 here!
! --------------------------------------------------------------------
        DO 20 I=1,3
          T   (I) = P(I) - POS(I,IAT)
          ILOW(I) = NINT( T(I) - ATMS(I,JATOM) + 0.4999D0 )
          IUP (I) = NINT( T(I) + ATMS(I,JATOM) - 0.4999D0 )
!         IF( ILOW(I) .GT. IUP(I) ) GOTO 10
   20   CONTINUE
!
!       << check all relevent atomic sphere displacements >>
        DO 30 JZ=ILOW(3),IUP(3)
          DO 40 JY=ILOW(2),IUP(2)
            DO 50 JX=ILOW(1),IUP(1)
!
!             << load x - R0 - R in Cartesian coordinates >>
              DO 60 J=1,3
                R(J) = (T(1)-JX)*BR2(1,J) + (T(2)-JY)*BR2(2,J) &
                     + (T(3)-JZ)*BR2(3,J)
   60         CONTINUE
              R2 = R(1)*R(1)+R(2)*R(2)+R(3)*R(3)
              IF( R2.LT.RMT2 ) THEN
!               << in muffin tin sphere IAT >>
                RR = SQRT(R2)
                ILAT(1) = JX
                ILAT(2) = JY
                ILAT(3) = JZ
                IR = MIN(1+INT(LOG(MAX(RR/stru%r0(JATOM),1.0D0))/DX(JATOM)) &
                        ,stru%Npt(JATOM))
                RETURN
              ENDIF
   50       CONTINUE
   40     CONTINUE
   30   CONTINUE
   10 CONTINUE
!     << in interstitial >>
      IAT = 0
      DO 70 J=1,3
         R(J) = P(1)*BR2(1,J) + P(2)*BR2(2,J) + P(3)*BR2(3,J)
   70 CONTINUE
!
      RETURN
      END


!!/---
!! Local Variables:
!! mode: f90
!! End:
!!\---
!!
!! Time-stamp: <2016-07-15 11:43:11 assman@faepop71.tu-graz.ac.at>
