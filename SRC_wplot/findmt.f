!!! wien2wannier/SRC_wplot/findmt.f
!!!
!!! $Id: findmt.f 385 2015-06-01 13:08:18Z assmann $

      SUBROUTINE FINDMT(P,ATMS,nato,NAT,IAT,ILAT,IR,R)
      use struct
      use radgrd
      use param
      use latt
      IMPLICIT REAL(R8) (A-H,O-Z)
      DIMENSION P(3),ATMS(3,NATO),ILAT(3),R(3)
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
      DO 10 IAT=1,NAT
        JATOM = IATNR(IAT)
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
                IR = MIN(1+INT(LOG(MAX(RR/RNOT(JATOM),1.0D0))/DX(JATOM)) &
                        ,JRI(JATOM))
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
!! Time-stamp: <2015-05-23 19:58:48 elias>
