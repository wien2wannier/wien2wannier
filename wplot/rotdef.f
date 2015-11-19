      SUBROUTINE ROTDEF(NAT,MULT,POS,IOP)
!     last changes: 08.11.00 ub (MVATOM introduced)
!                   12.11.00 ub (C:17/TR(I)CK17 feature with parameters)
!
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'param.inc'
      DIMENSION MULT(NAT), POS(3,Nat*48), IOP(Nat*48)
!
! << Input >>
! NAT      -- number of inequivalent basis atoms
! MULT(j)  -- number of equivalent atom within the j-th set of inequiv. atoms
! POS(:,i) -- primitive fractional coordinates of the i-th atom 
!
! from COMMON /SYM2/ -- symmetry operations
! IMAT,TAU -- symmetry operations (in primitive fractional coordinates)
! IORD     -- number of symmetry operations
! ----------------------------------------------------------------------------
! {Q|t} : y_i = Sum(j) Q_ij x_j + t_i  with  Q_ij = IMAT(j,i) and t_i = TAU(i)
! ----------------------------------------------------------------------------
!
! << Output >>
! POS(:,i) -- the proper symmetry generated position of the i-th atom
! IOP(i)   -- the symmetry operation which generates the i-th atom from
!             the first atom of the corresponding set of inequivalent atoms
!
! Caution: The symmetry operations selected here to generate the individual 
!          atoms must precisely co-incide with those selected in LAPW2 !
!
      COMMON /LATT/ VUC,BR1(3,3),BR2(3,3),BR3(3,3),BR4(3,3)
      COMMON /SYM2/ IMAT(3,3,NSYM),TAU(3,NSYM),IORD
      DIMENSION POS0(3),R(3),D(3)
!
      DATA TOL /1.0D-4/ , ONE,HALF /1.0D0, 0.5D0/
!
!     << find the generating symmetry operations >>
      write(*,*)"debug1rotdef"
      WRITE(6,2000)
      write(*,*)"debug2rotdef"
      INDEX=0
      DO 10 JATOM=1,NAT
        INDEX1 = INDEX+1
!       << store r_1.atom before updating it >>
        DO 15 J=1,3
          POS0(J) = POS(J,INDEX1)
   15   CONTINUE
        DO 20 M=1,MULT(JATOM)
          INDEX = INDEX+1
          DO 30 JOP=1,IORD
!
!           << find {Q|t}(r_1.atom) >>
            DO 40 I=1,3
              R(I) = IMAT(1,I,JOP)*POS0(1) &
                   + IMAT(2,I,JOP)*POS0(2) &
                   + IMAT(3,I,JOP)*POS0(3) + TAU(I,JOP)
   40       CONTINUE
!
!           << check the difference modulo lattice translations >>
            DO 50 I=1,3
              D(I) = ABS(MOD(ABS(POS(I,INDEX)-R(I))+HALF,ONE)-HALF)
   50       CONTINUE
            IF (D(1).LT.TOL.AND.D(2).LT.TOL.AND.D(3).LT.TOL) THEN
              IOP(INDEX) = JOP
!:17[
              IF(.NOT.MVATOM)THEN
                DO 55 I=1,3
   55           R(I) = POS(I,INDEX)
              ENDIF
!:17]
              DO 60 I=1,3
   60         POS(I,INDEX) = R(I)
!
!             << print updated position in Cartesian coordinates >>
              DO 70 J=1,3
                R(J) = POS(1,INDEX)*BR2(1,J) &
                     + POS(2,INDEX)*BR2(2,J) &
                     + POS(3,INDEX)*BR2(3,J)
   70         CONTINUE
              WRITE(6,2010) JATOM,IOP(INDEX),(R(J),J=1,3)
              GOTO 20
            END IF
   30     CONTINUE
!
!         << something is wrong here >>
          WRITE(6,1000) INDEX1,INDEX
          WRITE(6,1010) INDEX1,(POS(I,INDEX1),I=1,3)
          WRITE(6,1010) INDEX ,(POS(I,INDEX ),I=1,3)
          STOP 'ROTDEF'
!
 20     CONTINUE
 10   CONTINUE
!
      RETURN
!
 1000 FORMAT(///3X,'ERROR IN ROTDEF:'/ &
             'NO SYMMETRY OPERATION FOUND TO MAP ATOM',I3, &
             ' ONTO ATOM',I3)
 1010 FORMAT('ATOM',I3,' : POS = ',3F13.7)
!
 2000 FORMAT(/' SYMMETRY ADJUSTED POSITIONS OF THE BASIS ATOMS' &
             /' ----------------------------------------------' &
             /' atom  symm.     x [bohr]     y [bohr]     z [bohr]')
 2010 FORMAT(2(I5,1X),3F13.7)
!
      END
