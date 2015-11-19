!!! wien2wannier/SRC_wplot/grdgen.f
!!!
!!! $Id: grdgen.f 184 2014-02-04 20:34:24Z assmann $

      SUBROUTINE GRDGEN(MODE,NP,NPO)
      use grid
!
      use param
      use latt
      IMPLICIT REAL(R8) (A-H,O-Z)
      CHARACTER(4) MODE
      DIMENSION NP(3),NPO(0:3)
!
! Input Format
! ~~~~~~~~~~~~
! Record 1: (A4)
! MODE   --  'ANY '     if an arbitrary list of grid points is to be handled 
!            '<n>D <f>' if an n-dim. grid of grid points is to be handled
!                       n = 0, 1, 2, 3 and f = ' ', 'O'(ORTHO.), 'N'(ON-ORTHO.)
! 
! MODE == 'ANY '
! --------------
! Record 2: (*)
! NPG  --  total number of grid points
! NPO  --  echo output increment or zero for just echoing the first
!          and last grid point
!
! FOR EACH grid point
!   Record 3: from unit 7 (*) 
!   RGRID(1:3)  --  the grid points in Cartesian coordinates
! DONE
!
! MODE == '<n>D <f>'
! ------------------
! Record 2: (*)
! IX,IY,IZ,IDIV  --  the origin (IX/IDIV,IY/IDIV,IZ/IDIV) of the evaluation
!                    grid in conventional fractional coordinates
!
! FOR EACH dimension a = 1,...,n
!   Record 3: (*)
!   IX,IY,IZ,IDIV  --  the end-point (IX/IDIV,IY/IDIV,IZ/IDIV) of the a-axis
!                      in conventional fractional coordinates
! DONE
!
! Record 4: (*)
! NPX,NPY,...    --  number of grid points along each axis
! NPOX,NPOY,...  --  echo output increment along each axis or zero for just
!                    taking the first and last point on each axis
!
! If f = ' ' or 'O' the axes are checked for orthogonality !
!
! Output:
! ~~~~~~~
! MODE       -- input mode
!
! for MODE = 'ANY'
! NPO(0)     -- echo output increment
! 
! for MODE = '<n>D <f>' (n > 0)
! NP(a)      -- number of grid points along the a-axis
! NPO(a)     -- echo output increment along the a-axis
!
! in any case 
! NPG        -- total number of grid points
! RGRID(:,i) -- the i-th grid point in primitive fractional coordinates


      DIMENSION O(3),F(3,3),X(3,0:3),RAX(3),PHI(1:2,2:3)
      CHARACTER AXIS(3)
      DATA AXIS / 'x', 'y', 'z' /

      character(*), parameter :: fmt_3d = &
           "(' 3D-NET OF GRID POINTS'                     / &
           & ' ---------------------'                     / &
           & ' number of grid points for x, y, z: ', 3I4,   &
           & ' (total:', I7, ')'                          / &
           & ' output echo increment for x, y, z: ', 3I4  / &
           & ' number of echo points for x, y, z: ', 3I4,   &
           & ' (total:', I7, ')')"
      character(*), parameter :: fmt_2d = &
           "(' 2D-NET OF GRID POINTS'                     / &
           & ' ---------------------'                     / &
           & ' number of grid points for x, y: ', 2I5,      &
           & ' (total:',I7, ')'                           / &
           & ' output echo increment for x, y: ', 2I5     / &
           & ' number of echo points for x, y: ', 2I5,      &
           & ' (total:',I7, ')')"
      character(*), parameter :: fmt_1d = &
           "(' 1D-NET OF GRID POINTS'                     / &
           & ' ---------------------'                     / &
           & ' number of grid points: ',I7                / &
           & ' output echo increment: ',I7                / &
           & ' number of echo points: ',I7)"
      character(*), parameter :: fmt_0d = &
           "(' 0D-NET OF GRID POINTS'                     / &
           & ' ---------------------')"

      character(512) :: fmt_grid

!
!     <<  in the plotting grid >>
      READ(unit_in,1000) MODE
      IF(MODE(1:3).EQ.'ANY')THEN
!       << arbitray list >>
        READ(unit_in,*) NPG,NPO(0)

        ngdim=npg
        allocate (RGRID(3,NGDIM),IREG(NGDIM),ILAT(3,NGDIM),IRI(NGDIM))
!        IF(NPG.GT.NGDIM) STOP 'NGDIM too small (1)'
        IF(NPO(0).EQ.0) NPO(0)=MAX(NPG-1,1)
        WRITE(unit_out,2000) NPG,NPO(0),1+(NPG-1)/NPO(0)
!       << read in grid points >>
        DO 10 IG=1,NPG
          READ(unit_grid,*) (O(J),J=1,3)
!         << transform to primitive fractional coordinates >>
          DO 20 I=1,3
            RGRID(I,IG) = BR4(I,1)*O(1)+BR4(I,2)*O(2)+BR4(I,3)*O(3)
   20     CONTINUE
   10   CONTINUE
!       << write grid info on data file >>
        WRITE(unit_psink,3000) MODE,NPG
      ELSE
!       << the origin >>
        READ(unit_in,*) (O(I),I=1,3),IDV

        O = O / IDV
        NDIM = 0
        fmt_grid = fmt_0d
        IF(MODE(1:2).NE.'0D')THEN
!         << x-end >>
          READ(unit_in,*) (F(I,1),I=1,3),IDV

          F(:,1) = F(:,1) / IDV - O(:)
          NDIM = 1
          fmt_grid = fmt_1d
          IF(MODE(1:2).NE.'1D')THEN
!           << y-end >>
            READ(unit_in,*) (F(I,2),I=1,3),IDV
            F(:,2) = F(:,2) / IDV - O(:)
            NDIM = 2
            fmt_grid = fmt_2d
            IF(MODE(1:2).NE.'2D')THEN
!             << z-end >>
              READ(unit_in,*) (F(I,3),I=1,3),IDV
              F(:,3) = F(:,3) / IDV - O(:)
              NDIM = 3
              fmt_grid = fmt_3d
            ENDIF
          ENDIF
        ENDIF

!       << grid size and echo increments >>
        IF(NDIM.GT.0) READ(unit_in,*) (NP(K),K=1,NDIM),(NPO(K),K=1,NDIM)

        close(unit_in)

        NPG = 1
        NPGO = 1
        DO 70 K=1,NDIM
          NPG = NPG * NP(K)
          IF(NPO(K).EQ.0) NPO(K)=MAX(NP(K)-1,1)
          NPGO = NPGO * ( 1+(NP(K)-1)/NPO(K) )
   70   CONTINUE
        DO 80 K=NDIM+1,3
          NP (K) = 1
          NPO(K) = 1
   80   CONTINUE
        IF( NDIM.GT.1 )THEN
          WRITE(unit_out,fmt_grid) (NP(K),K=1,NDIM),NPG,(NPO(K),K=1,NDIM), &
                         (1+(NP(K)-1)/NPO(K),K=1,NDIM),NPGO
        ELSE IF( NDIM.EQ.1 )THEN
          WRITE(unit_out,fmt_grid) NP(1),NPO(1),1+(NP(1)-1)/NPO(1)
        ELSE
          WRITE(unit_out,fmt_grid)
        ENDIF
!
!       << transform origin and axes into primitive fractional coordinates >>
!       << a) transform into Cartesian coordinates >>
!       <<    and check axes for orthogonality     >>
        IF(MODE(4:4).EQ.' ')MODE(4:4)='O'
        WRITE(unit_out,2050)
        X(:,0) = O(1)*BR1(1,:)+O(2)*BR1(2,:)+O(3)*BR1(3,:)
        WRITE(unit_out,2060) (O(I),I=1,3),(X(J,0),J=1,3)
        DO 100 K=1,NDIM
          DO 110 J=1,3
            X(J,K) = F(1,K)*BR1(1,J)+F(2,K)*BR1(2,J)+F(3,K)*BR1(3,J)
  110     CONTINUE
          WRITE(unit_out,2070) AXIS(K),(F(I,K),I=1,3),(X(J,K),J=1,3)
  100   CONTINUE
        WRITE(unit_out,2080)
        DO 115 K=1,NDIM
          RAX(K) = SQRT( X(1,K)*X(1,K)+X(2,K)*X(2,K)+X(3,K)*X(3,K) )
          WRITE(unit_out,2090) AXIS(K),RAX(K)
  115   CONTINUE
        IF(NDIM.GT.1)WRITE(unit_out,2100)
        DO I=1,NDIM-1
          DO J=I+1,NDIM
            RIJ = X(1,I)*X(1,J) + X(2,I)*X(2,J) + X(3,I)*X(3,J)
            COSPHI = RIJ/(RAX(I)*RAX(J))
            PHI(I,J) = ACOS(COSPHI)*57.295779513082321D0
            WRITE(unit_out,2110) AXIS(I),AXIS(J),COSPHI,PHI(I,J)
            IF(MODE(4:4).EQ.'O'.AND.ABS(COSPHI).GT.0.001) &
               STOP 'NON-ORTHOGONAL AXES'
         END DO
      END DO
!       << b) transform into primitive fractional coordinates >>
      O(:) = BR4(:,1)*X(1,0)+BR4(:,2)*X(2,0)+BR4(:,3)*X(3,0)
      DO K=1,NDIM
         DO I=1,3
            F(I,K) = BR4(I,1)*X(1,K)+BR4(I,2)*X(2,K)+BR4(I,3)*X(3,K)
         END DO
      END DO
!
!       << generate the evaluation grid >>
        ngdim=npg
        allocate (RGRID(3,NGDIM),IREG(NGDIM),ILAT(3,NGDIM),IRI(NGDIM))
!        IF(NPG.GT.NGDIM) STOP 'NGDIM too small (2)'
        DO I=1,3
           DO K=1,NDIM
              F(I,K) = F(I,K) / MAX(NP(K)-1,1)
           END DO
           DO K=NDIM+1,3
              F(I,K) = 0.0D0
           END DO
        END DO
!       << use gnuplot order for data generation, i.e.  >>
!       << POS = IZ + (IY-1)*NP(3) + (IX-1)*NP(3)*NP(2) >>
        IG = 0
        DO 190 IX=0,NP(1)-1
          DO 200 IY=0,NP(2)-1
            DO 210 IZ=0,NP(3)-1
              IG = IG + 1
              DO 220 I=1,3
                RGRID(I,IG) = O(I) + IX*F(I,1)+IY*F(I,2)+IZ*F(I,3)
  220         CONTINUE
  210       CONTINUE
  200     CONTINUE
  190   CONTINUE
!       << write grid info on data file >>
        MODE(4:4)=' '
        IF(NDIM.EQ.0)THEN
          WRITE(unit_psink,3010)MODE
        ELSE
          WRITE(unit_psink,3020)MODE,(J,J=1,NDIM-1)
          DO 127 I=1,NDIM
            WRITE(unit_psink,3030)NP(I),RAX(I),(PHI(J,I),J=1,I-1)
  127     CONTINUE
          IF(NDIM.EQ.2)WRITE(unit_psink,3040)
          IF(NDIM.EQ.3)WRITE(unit_psink,3045)
        ENDIF
      ENDIF
    RETURN
!
 1000 FORMAT(A4)
!
 2000 FORMAT( ' ARBITRARY LIST OF GRID POINTS' &
             /' -----------------------------' &
             /' number of grid points: ',I7 &
             /' output echo increment: ',I7 &
             /' number of echo points: ',I7)
 2050 FORMAT(/' PLOTTING AREA' &
             /' -------------' &
             /' x = Sum(j=1,3) f_i a_i  with  f_i in', &
              ' conventional fractional coordinates' &
            //'           f_1      f_2      f_3   ', &
              '        x [bohr]     y [bohr]     z [bohr]')
 2060 FORMAT(' origin',1X,3F9.5,3X,3F13.7)
 2070 FORMAT(' ',A1,'-axis',1X,3F9.5,3X,3F13.7)
 2080 FORMAT(/'        length [bohr]')
 2090 FORMAT(' abs(',A1,') ',F13.7)
 2100 FORMAT(/'          cos(phi)   phi [degree]')
 2110 FORMAT(' <)(',A1,',',A1,') ',F10.7,3X,F8.3)
!
 3000 FORMAT(A4,'  NP =',I8,/, &
             'order according to the grid points ', &
             'provided in the input file <case>.grid')
 3010 FORMAT(A4)
 3020 FORMAT(A4,'  NP         abs(X)': &
             '   ang(X',I1,',X)':'   ang(X',I1,',X)')
 3030 FORMAT(I8,2X,F13.7,2(2X,F10.5))
 3040 FORMAT('order: ((psi(ix,iy),iy=1,ny),ix=1,nx)')
 3045 FORMAT('order: (((psi(ix,iy,iz),iz=1,nz),iy=1,ny),ix=1,nx)')
      END

!!/---
!! Local Variables:
!! mode: f90
!! End:
!!\---
