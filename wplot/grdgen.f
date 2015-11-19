      SUBROUTINE GRDGEN(MODE,NP,NPO)
      use grid
!     last changes: 19.07.00 ba (radfu2 stuff)
!                   15.08.00 ub (merging and bug-fixing)
!                   31.08.00 ub (update COMMON /GRID/)
!                   02.11.00 ub (optimizing the user interface)
!                   13.11.00 ub (restore COMMON /GRID/)
!
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'param.inc'
      CHARACTER MODE*4
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
!
!      COMMON /GRID  / RGRID(3,NGDIM),NPG, &
!                      IREG(NGDIM),ILAT(3,NGDIM),IRI(NGDIM)
      COMMON /LATT/ VUC,BR1(3,3),BR2(3,3),BR3(3,3),BR4(3,3)
      DIMENSION O(3),F(3,3),X(3,0:3),RAX(3),PHI(1:2,2:3)
      CHARACTER*1 AXIS(3)
      DATA AXIS / 'x', 'y', 'z' /
!
!     <<  in the plotting grid >>
      READ(5,1000) MODE
      IF(MODE(1:3).EQ.'ANY')THEN
!       << arbitray list >>
        READ(5,*) NPG,NPO(0)
        !write(*,*)"debuggrid",NPG
        ngdim=npg
        allocate (RGRID(3,NGDIM),IREG(NGDIM),ILAT(3,NGDIM),IRI(NGDIM))
!        IF(NPG.GT.NGDIM) STOP 'NGDIM too small (1)'
        IF(NPO(0).EQ.0) NPO(0)=MAX(NPG-1,1)
        WRITE(6,2000) NPG,NPO(0),1+(NPG-1)/NPO(0)
!       << read in grid points >>
        DO 10 IG=1,NPG
          READ(7,*) (O(J),J=1,3)
!         << transform to primitive fractional coordinates >>
          DO 20 I=1,3
            RGRID(I,IG) = BR4(I,1)*O(1)+BR4(I,2)*O(2)+BR4(I,3)*O(3)
   20     CONTINUE
   10   CONTINUE
!       << write grid info on data file >>
        WRITE(21,3000) MODE,NPG
      ELSE
!       << the origin >>
        READ(5,*) (O(I),I=1,3),IDV
        !write(*,*)"debug orig",O
        DO 30 I=1,3
   30   O(I) = O(I) / DBLE(IDV)
        NDIM = 0
        assign 2040 to IFORM
        IF(MODE(1:2).NE.'0D')THEN
!         << x-end >>
          READ(5,*) (F(I,1),I=1,3),IDV
          !write(*,*)"debug vec1",F(:,1)
          DO 40 I=1,3
   40     F(I,1) = F(I,1) / DBLE(IDV) - O(I)
          NDIM = 1
          assign 2030 to IFORM
          IF(MODE(1:2).NE.'1D')THEN
!           << y-end >>
            READ(5,*) (F(I,2),I=1,3),IDV
            DO 50 I=1,3
   50       F(I,2) = F(I,2) / DBLE(IDV) - O(I)
            NDIM = 2
            assign 2020 to IFORM
            IF(MODE(1:2).NE.'2D')THEN
!             << z-end >>
              READ(5,*) (F(I,3),I=1,3),IDV
              DO 60 I=1,3
   60         F(I,3) = F(I,3) / DBLE(IDV) - O(I)
              NDIM = 3
              assign 2010 to IFORM
            ENDIF
          ENDIF
        ENDIF
        !write(*,*)"debug",F
!       << grid size and echo increments >>
        IF(NDIM.GT.0) READ(5,*) (NP(K),K=1,NDIM),(NPO(K),K=1,NDIM)
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
          WRITE(6,IFORM) (NP(K),K=1,NDIM),NPG,(NPO(K),K=1,NDIM), &
                         (1+(NP(K)-1)/NPO(K),K=1,NDIM),NPGO
        ELSE IF( NDIM.EQ.1 )THEN
          WRITE(6,IFORM) NP(1),NPO(1),1+(NP(1)-1)/NPO(1)
        ELSE
          WRITE(6,IFORM)
        ENDIF
!
!       << transform origin and axes into primitive fractional coordinates >>
!       << a) transform into Cartesian coordinates >>
!       <<    and check axes for orthogonality     >>
        IF(MODE(4:4).EQ.' ')MODE(4:4)='O'
        WRITE(6,2050)
        DO 90 J=1,3
   90   X(J,0) = O(1)*BR1(1,J)+O(2)*BR1(2,J)+O(3)*BR1(3,J)
        WRITE(6,2060) (O(I),I=1,3),(X(J,0),J=1,3)
        DO 100 K=1,NDIM
          DO 110 J=1,3
            X(J,K) = F(1,K)*BR1(1,J)+F(2,K)*BR1(2,J)+F(3,K)*BR1(3,J)
  110     CONTINUE
          WRITE(6,2070) AXIS(K),(F(I,K),I=1,3),(X(J,K),J=1,3)
  100   CONTINUE
        WRITE(6,2080)
        DO 115 K=1,NDIM
          RAX(K) = SQRT( X(1,K)*X(1,K)+X(2,K)*X(2,K)+X(3,K)*X(3,K) )
          WRITE(6,2090) AXIS(K),RAX(K)
  115   CONTINUE
        IF(NDIM.GT.1)WRITE(6,2100)
        DO 120 I=1,NDIM-1
          DO 125 J=I+1,NDIM
            RIJ = X(1,I)*X(1,J) + X(2,I)*X(2,J) + X(3,I)*X(3,J)
            COSPHI = RIJ/(RAX(I)*RAX(J))
            PHI(I,J) = ACOS(COSPHI)*57.295779513082321D0
            WRITE(6,2110) AXIS(I),AXIS(J),COSPHI,PHI(I,J)
            IF(MODE(4:4).EQ.'O'.AND.ABS(COSPHI).GT.0.001) &
               STOP 'NON-ORTHOGONAL AXES'
  125     CONTINUE
  120   CONTINUE
!       << b) transform into primitive fractional coordinates >>
        DO 130 I=1,3
  130   O(I) = BR4(I,1)*X(1,0)+BR4(I,2)*X(2,0)+BR4(I,3)*X(3,0)
        DO 140 K=1,NDIM
          DO 150 I=1,3
  150     F(I,K) = BR4(I,1)*X(1,K)+BR4(I,2)*X(2,K)+BR4(I,3)*X(3,K)
  140   CONTINUE
!
!       << generate the evaluation grid >>
        ngdim=npg
        allocate (RGRID(3,NGDIM),IREG(NGDIM),ILAT(3,NGDIM),IRI(NGDIM))
!        IF(NPG.GT.NGDIM) STOP 'NGDIM too small (2)'
        DO 160 I=1,3
          DO 170 K=1,NDIM
  170     F(I,K) = F(I,K) / MAX(NP(K)-1,1)
          DO 180 K=NDIM+1,3
  180     F(I,K) = 0.0D0
  160   CONTINUE
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
          WRITE(21,3010)MODE
        ELSE
          WRITE(21,3020)MODE,(J,J=1,NDIM-1)
          DO 127 I=1,NDIM
            WRITE(21,3030)NP(I),RAX(I),(PHI(J,I),J=1,I-1)
  127     CONTINUE
          IF(NDIM.EQ.2)WRITE(21,3040)
          IF(NDIM.EQ.3)WRITE(21,3045)
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
 2010 FORMAT( ' 3D-NET OF GRID POINTS' &
             /' ---------------------' &
             /' number of grid points for x, y, z: ',3I4,' (total:',I7, &
          ')'/' output echo increment for x, y, z: ',3I4 &
             /' number of echo points for x, y, z: ',3I4,' (total:',I7, &
          ')')
 2020 FORMAT( ' 2D-NET OF GRID POINTS' &
             /' ---------------------' &
             /' number of grid points for x, y: ',2I5,' (total:',I7, &
          ')'/' output echo increment for x, y: ',2I5 &
             /' number of echo points for x, y: ',2I5,' (total:',I7, &
          ')')
 2030 FORMAT( ' 1D-NET OF GRID POINTS' &
             /' ---------------------' &
             /' number of grid points: ',I7 &
             /' output echo increment: ',I7 &
             /' number of echo points: ',I7)
 2040 FORMAT( ' 0D-NET OF GRID POINTS' &
             /' ---------------------')
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
