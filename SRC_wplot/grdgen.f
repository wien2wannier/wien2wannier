!!! wien2wannier/SRC_wplot/grdgen.f
!!!
!!!    Generates / reads real-space grid for ‘wplot’.

SUBROUTINE GRDGEN(MODE,NP)

  use const, only: BUFSZ
  use grid,  only: R8, NPG, Rgrid, ireg, ilat, IRI
  use param, only: unit_in, unit_out, unit_grid, unit_psink
  use latt,  only: BR1, BR4
  use clio,  only: croak
  use util,  only: string

  implicit none


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
!
! If f = ' ' or 'O' the axes are checked for orthogonality !
!
! Output:
! ~~~~~~~
! MODE       -- input mode
!
! for MODE = 'ANY'
! 
! for MODE = '<n>D <f>' (n > 0)
! NP(a)      -- number of grid points along the a-axis
!
! in any case 
! NPG        -- total number of grid points
! RGRID(:,i) -- the i-th grid point in primitive fractional coordinates


  character(4),intent(out) :: mode
  integer,     intent(out) :: NP(3)

  real(R8) :: O(3), F(3,3), X(3,0:3), RAX(3), PHI(1:2,2:3), Rij, cosphi
  integer  :: Ngdim, ig, i, j, idv, k, Ndim, ix, iy, iz

  character, parameter :: axis(3) = (/ 'x', 'y', 'z' /)

!!!   Format parameters (cannot be ‘parameter’ due to different
!!!   lengths)
  character(len=BUFSZ) :: fmt_grid(-1:3)

  fmt_grid(-1) = &
       "(' ARBITRARY LIST OF GRID POINTS' / &
       & ' -----------------------------' / &
       & ' number of grid points: ',I7)"
  fmt_grid(0) = &
       "(' 0D: 1 POINT'    / &
       & ' -----------')"
  fmt_grid(1) = &
       "(' 1D-NET OF GRID POINTS'         / &
       & ' ---------------------'         / &
       & ' number of grid points: ',I7)"

  fmt_grid(2) = &
       "(' 2D-NET OF GRID POINTS'                                      / &
       & ' ---------------------'                                      / &
       & ' number of grid points for x, y: ',2I5,' (total:',I7, ')')"
  fmt_grid(3) = &
       "(' 3D-NET OF GRID POINTS'                                         / &
       & ' ---------------------'                                         / &
       & ' number of grid points for x, y, z: ',3I4,' (total:', I7, ')')"


!!!   <<  in the plotting grid >>
  read(unit_in, '(A4)') MODE
  if(MODE(1:3) == 'ANY') then
     Ndim = -1

     ! << arbitray list >>
     read(unit_in,*) NPG

     ngdim=npg
     allocate (RGRID(3,NGDIM),IREG(NGDIM),ILAT(3,NGDIM),IRI(NGDIM))

     read_grid: do ig=1,NPG
        select case (MODE(4:4))
        case ('F')
           read(unit_grid, *) Rgrid(:, ig), IDV
           Rgrid(:, ig) = Rgrid(:, ig) / IDV

        case (' ', 'C')
           read(unit_grid,*) Rgrid(:, ig)
           ! << transform to primitive fractional coordinates >>
           Rgrid(:, ig) = matmul(BR4(:,:), Rgrid(:, ig))

        case default
           call croak('unknown FLAG: '//mode(4:4))
        end select
     end DO read_grid

     ! << write grid info on data file >>
     write(unit_psink,                                &
          "(A4,'  NP =',I8,/,                         &
          & 'order according to the grid points ',      &
          & 'provided in the input file <case>.grid')") &
          MODE,NPG
  else
     ! << the origin >>
     READ(unit_in,*) O, IDV
     O = O / IDV

     NDIM = 0
     IF(MODE(1:2).NE.'0D')THEN
        ! << x-end >>
        READ(unit_in,*) (F(I,1),I=1,3),IDV
        F(:,1) = F(:,1) / IDV - O

        NDIM = 1
        IF(MODE(1:2).NE.'1D')THEN
           ! << y-end >>
           READ(unit_in,*) (F(I,2),I=1,3),IDV
           F(:,2) = F(:,2) / IDV - O

           NDIM = 2
           IF(MODE(1:2).NE.'2D')THEN
              ! << z-end >>
              READ(unit_in,*) (F(I,3),I=1,3),IDV
              F(:,3) = F(:,3) / IDV - O

              NDIM = 3
           ENDIF
        ENDIF
     ENDIF

     ! << grid size >>
     IF(NDIM.GT.0) READ(unit_in,*) NP(1:NDIM)

     NPG = product(NP(1:Ndim))

     NP(NDIM+1 : 3) = 1

     
!!! Write grid info
     select case (Ndim)
     case (-1)
        write(unit_out, fmt_grid(Ndim)) NPG
     case (0)
        write(unit_out, fmt_grid(Ndim))
     case (1)
        write(unit_out, fmt_grid(Ndim)) NP(1)
     case (2, 3)
        write(unit_out, fmt_grid(Ndim)) NP(1:NDIM), NPG
     case default
        call croak("What is this, string theory?  Ndim="//&
             &     trim(string(Ndim)))
     end select

     ! << transform origin and axes into primitive fractional coordinates >>
     ! << a) transform into Cartesian coordinates >>
     ! <<    and check axes for orthogonality     >>
     select case (mode(4:4))
     case ('N', 'O')
        continue
     case (' ')
        mode(4:4) = 'O'
     case default
        call croak('unknown FLAG: '//mode(4:4))
     end select

     write(unit_out, "(/' PLOTTING AREA'                         &
          &            /' -------------'                         &
          &            /' x = Sum(j=1,3) f_i a_i  with  f_i in', &
          &             ' conventional fractional coordinates'   &
          &           //'           f_1      f_2      f_3   ',   &
          &             '        x [bohr]     y [bohr]     z [bohr]')")

     X(:,0) = O(1)*BR1(1,:)+O(2)*BR1(2,:)+O(3)*BR1(3,:)

     write(unit_out, "(' origin', 1X, 3F9.5, 3X, 3F13.7)") O, X(:,0)

     DO K=1,NDIM
        DO J=1,3
           X(J,K) = F(1,K)*BR1(1,J)+F(2,K)*BR1(2,J)+F(3,K)*BR1(3,J)
        end DO

        write(unit_out, "(' ',A1,'-axis',1X,3F9.5,3X,3F13.7)") &
             AXIS(K),F(:,K), X(:,K)
     end DO

     write(unit_out, "(/'        length [bohr]')")

     do K=1,NDIM
        RAX(K) = SQRT( X(1,K)*X(1,K)+X(2,K)*X(2,K)+X(3,K)*X(3,K) )
        write(unit_out, "(' abs(',A1,') ',F13.7)") AXIS(K), RAX(K)
     end do

     IF(NDIM>1) write(unit_out, "(/'          cos(phi)   phi [degree]')")
     do I=1,NDIM-1
        do J=I+1,NDIM
           RIJ = X(1,I)*X(1,J) + X(2,I)*X(2,J) + X(3,I)*X(3,J)
           COSPHI = RIJ/(RAX(I)*RAX(J))
           PHI(I,J) = ACOS(COSPHI)*57.295779513082321D0

           write(unit_out, "(' <)(',A1,',',A1,') ',F10.7,3X,F8.3)") &
                AXIS(I), AXIS(J), COSPHI, PHI(I,J)

           IF(MODE(4:4) == 'O' .AND. ABS(COSPHI)>0.001) &
                STOP 'NON-ORTHOGONAL AXES'
        end do
     end do
     
     ! << b) transform into primitive fractional coordinates >>
     O(:) = BR4(:,1)*X(1,0)+BR4(:,2)*X(2,0)+BR4(:,3)*X(3,0)
     DO K=1,NDIM
        DO I=1,3
           F(I,K) = BR4(I,1)*X(1,K)+BR4(I,2)*X(2,K)+BR4(I,3)*X(3,K)
        END DO
     END DO

     ! << generate the evaluation grid >>
     ngdim=npg
     allocate (RGRID(3,NGDIM),IREG(NGDIM),ILAT(3,NGDIM),IRI(NGDIM))
     DO I=1,3
        DO K=1,NDIM
           F(I,K) = F(I,K) / MAX(NP(K)-1,1)
        END DO
        DO K=NDIM+1,3
           F(I,K) = 0.0D0
        END DO
     END DO
     ! << use gnuplot order for data generation, i.e.  >>
     ! << POS = IZ + (IY-1)*NP(3) + (IX-1)*NP(3)*NP(2) >>
     ig = 0
     do ix=0,NP(1)-1
        do iy=0,NP(2)-1
           do iz=0,NP(3)-1
              ig = ig + 1
              Rgrid(:,ig) = O(:) + ix*F(:,1) + iy*F(:,2) + iz*F(:,3)
           end do
        end do
     end do
     ! << write grid info on data file >>
     MODE(4:4)=' '
     if(NDIM.eq.0)then
        WRITE(unit_psink, "(A4)") MODE
     else
        write(unit_psink,                                 &
             "(A4, '  NP         abs(X)':                 &
             & '   ang(X',I1,',X)':'   ang(X',I1,',X)')") &
             MODE, (J,J=1,NDIM-1)
        do I=1,NDIM
           WRITE(unit_psink, "(I8,2X,F13.7,2(2X,F10.5))") &
                NP(I),RAX(I),(PHI(J,I),J=1,I-1)
        end do
        if(NDIM == 2) &
             write(unit_psink, "('order: ((psi(ix,iy),iy=1,ny),ix=1,nx)')")
        if(NDIM == 3) &
             write(unit_psink, &
             "('order: (((psi(ix,iy,iz),iz=1,nz),iy=1,ny),ix=1,nx)')")
     endif
  endif
end SUBROUTINE GRDGEN

!!/---
!! Local Variables:
!! mode: f90
!! End:
!!\---
!!
!! Time-stamp: <2015-08-13 16:15:10 elias@hupuntu>
