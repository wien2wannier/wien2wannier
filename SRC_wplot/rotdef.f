!!! wien2wannier/SRC_wplot/rotdef.f

subroutine rotdef(stru, IOP)
  use const,     only: DPk
  use structmod, only: struct_t
  use clio,      only: croak
  use param,     only: unit_out, mvatom
  use latt,      only: br2
  use struct,    only: pos
  use sym2,      only: imat, trans, iord

  implicit none

  type(struct_t), intent(in)  :: stru
  integer,        intent(out) :: IOP(stru%Nneq*48)

  real(DPk), parameter :: half = 0.5_DPk, tol = 1e-4_DPk

  integer   :: index, index1, i, j, jatom, jOp, m
  real(DPk) :: pos0(3), R(3), D(3)

! << Input >>
! stru     -- struct
!
! from MODULE SYM2 -- symmetry operations
! IMAT,TRANS -- symmetry operations (in primitive fractional coordinates)
! IORD     -- number of symmetry operations
! ----------------------------------------------------------------------------
! {Q|t} : y_i = Sum(j) Q_ij x_j + t_i  with  Q_ij = IMAT(j,i) and t_i = TRANS(i)
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
!     << find the generating symmetry operations >>
      WRITE(unit_out,2000)

      INDEX=0
      DO 10 JATOM=1,STRU%NNEQ
        INDEX1 = INDEX+1
!       << store r_1.atom before updating it >>
        DO 15 J=1,3
          POS0(J) = STRU%POS(J,INDEX1)
   15   CONTINUE
        do 20 M=1,stru%mult(jatom)
          INDEX = INDEX+1
          DO 30 JOP=1,IORD
!
!           << find {Q|t}(r_1.atom) >>
            DO 40 I=1,3
              R(I) = IMAT(1,I,JOP)*POS0(1) &
                   + IMAT(2,I,JOP)*POS0(2) &
                   + IMAT(3,I,JOP)*POS0(3) + TRANS(I,JOP)
   40       CONTINUE
!
!           << check the difference modulo lattice translations >>
            DO 50 I=1,3
               D(I) = abs(mod(abs(stru%pos(i,index)-R(I)) + HALF, 1._DPk) &
                    - HALF)
   50       CONTINUE
            IF (D(1).LT.TOL.AND.D(2).LT.TOL.AND.D(3).LT.TOL) THEN
              IOP(INDEX) = JOP
!:17[
              IF(.NOT.MVATOM)THEN
                DO I=1,3
                   R(I) = STRU%POS(I,INDEX)
                END DO
              ENDIF
!:17]
              DO I=1,3
                 POS(I,INDEX) = R(I)
              END DO
!
!             << print updated position in Cartesian coordinates >>
              DO 70 J=1,3
                R(J) = POS(1,INDEX)*BR2(1,J) &
                     + POS(2,INDEX)*BR2(2,J) &
                     + POS(3,INDEX)*BR2(3,J)
   70         CONTINUE
              WRITE(unit_out,2010) JATOM,IOP(INDEX),(R(J),J=1,3)
              GOTO 20
            END IF
   30     CONTINUE
!
!         << something is wrong here >>
          WRITE(unit_out,1000) INDEX1,INDEX
          WRITE(unit_out,1010) INDEX1,(POS(I,INDEX1),I=1,3)
          WRITE(unit_out,1010) INDEX ,(POS(I,INDEX ),I=1,3)
          call croak('error in ROTDEF')
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


!!/---
!! Local Variables:
!! mode: f90
!! End:
!!\---
!!
!! Time-stamp: <2016-07-15 11:41:39 assman@faepop71.tu-graz.ac.at>
