!!! wien2wannier/SRC_wplot/latgen.f

subroutine LATGEN(stru)
  use const,     only: DPk, TAU, SQ3
  use latt,      only: br1, br2, br3, br4
  use param,     only: unit_out
  use structmod, only: struct_t
  use clio,      only: croak

  implicit none

  type(struct_t), intent(in) :: stru

  real(DPk)    :: alpha, beta, gamma, cosPhi, phi
  integer      :: i, j

! << Output (to module LATT) >>
! BR1(i,:) -- the real space lattice vectors a_i of the conventional unit cell
! BR2(i,:) -- the real space lattice vectors a_i of the primitive unit cell
! BR3(i,:) -- the reciprocal lattice vectors b_i of the conventional unit cell
! BR4(i,:) -- the reciprocal lattice vectors b_i of the primitive unit cell
!
! Here, the reciprocal lattice vectors are evaluated without a factor of 2pi !
!
! Caution: The lattice vectors setup here must precisely co-incide with
!          those used within LAPW2 !
!
      ALPHA = stru%alpha(1) * TAU/360
      BETA  = stru%alpha(2) * TAU/360
      GAMMA = stru%alpha(3) * TAU/360

      BR1=0; BR2=0

      IF(STRU%LATTIC(1:1).EQ.'P')THEN
! -------------------------------------------------------------------------
!       << primitive lattice : P a b c alp bet gam >>
!
!       a_1 = a * (sin(bet)*sin(phi),sin(bet)*cos(phi),cos(bet))
!       a_2 = b * (        0        ,      sin(alp)   ,cos(alp))
!       a_3 = c * (        0        ,        0        ,  1     )
!       with
!       cos(phi) := ( cos(gam) - cos(alp)*cos(bet) ) / sin(alp)*sin(bet)
!
!       triclinic, monoclinic, orthorhombic, tetragonal, cubic
!
!       b_1 = ( 1/sin(bet)*sin(phi)                    ,     0     ,0) / a
!       b_2 = (-1/sin(alp)*tan(phi)                    , 1/sin(alp),0) / b
!       b_3 = ( 1/tan(alp)*tan(phi)-1/tan(bet)*sin(phi),-1/tan(alp),1) / c
! -------------------------------------------------------------------------
        COSPHI=(COS(GAMMA)-COS(ALPHA)*COS(BETA))/SIN(ALPHA)/SIN(BETA)
        PHI=ACOS(COSPHI)
!       << primitive unit cell >>
        BR2(1,1)=STRU%A(1)*SIN(BETA)*SIN(PHI)
        BR2(1,2)=STRU%A(1)*SIN(BETA)*COS(PHI)
        BR2(1,3)=STRU%A(1)*COS(BETA)
        BR2(2,2)=STRU%A(2)*SIN(ALPHA)
        BR2(2,3)=STRU%A(2)*COS(ALPHA)
        BR2(3,3)=STRU%A(3)
!       << conventional unit cell >>
        BR1 = BR2
      ELSE IF(STRU%LATTIC(1:1).EQ.'H') THEN
! -------------------------------------------------------------------------
!       << hexagonal lattice : H a * c * * * >>
!
!       a_1 = a * (sqrt(3)/2,-1/2,0)
!       a_2 = a * (     0   ,  1 ,0)
!       a_3 = c * (     0   ,  0 ,1)
!
!       this setting corresponds to the P a a c 90 90 120 setting
!
!       b_1 = ( 1/sqrt(3),  -1 ,   1 ) / a
!       b_2 = ( 1/sqrt(3),   1 ,   1 ) / a
!       b_3 = (-2/sqrt(3),   0 ,   1 ) / c
! -------------------------------------------------------------------------
!       << primitive unit cell >>
        BR2(1,1)= STRU%A(1)*SQ3/2
        BR2(1,2)=-STRU%A(1)/2
        BR2(2,2)= STRU%A(1)
        BR2(3,3)= STRU%A(3)
!       << conventional unit cell >>
        BR1 = BR2
      ELSE IF(STRU%LATTIC(1:1).EQ.'T'.OR.STRU%LATTIC(1:1).EQ.'R') THEN
! -------------------------------------------------------------------------
!       << rhombohedral or trigonal lattice : R a * c * * * >>
!
!       a_1 = ( a/(2*sqrt(3)),-a/2,c/3)
!       a_2 = ( a/(2*sqrt(3)), a/2,c/3)
!       a_3 = (-a/   sqrt(3) ,  0 ,c/3)
!
!       Note: Although the trigonal lattice is treated as a primitive
!             lattice the lattice parameter correspond to the surrounding
!             (non-primitive) hexagonal unit cell.
!
!       b_1 = ( 1/(a*sqrt(3)), -1/b , 1/c )
!       b_2 = ( 1/(a*sqrt(3)),  1/b , 1/c )
!       b_3 = (-2/(a*sqrt(3)),   0  , 1/c )
! -------------------------------------------------------------------------
!       << primitive unit cell >>
        BR2(1,1)= STRU%A(1)/2/SQ3
        BR2(1,2)=-STRU%A(1)/2
        BR2(1,3)= STRU%A(3)/3
        BR2(2,1)= STRU%A(1)/2/SQ3
        BR2(2,2)= STRU%A(1)/2
        BR2(2,3)= STRU%A(3)/3
        BR2(3,1)=-STRU%A(1)/SQ3
        BR2(3,3)= STRU%A(3)/3
!       << conventional unit cell >>
        BR1 = BR2
      ELSE IF(STRU%LATTIC(1:1).EQ.'F') THEN
! -------------------------------------------------------------------------
!       << face-centered lattice : F a b c * * *
!
!       a_1 = ( 0 ,b/2,c/2)
!       a_2 = (a/2, 0 ,c/2)
!       a_3 = (a/2,b/2, 0 )
!
!       orthorhombic, cubic
!
!       b_1 = ( -1/a ,  1/b ,  1/c )
!       b_2 = (  1/a , -1/b ,  1/c )
!       b_3 = (  1/a ,  1/b , -1/c )
! -------------------------------------------------------------------------
!       << primitive unit cell >>
        BR2(1,2)=STRU%A(2)/2
        BR2(1,3)=STRU%A(3)/2
        BR2(2,1)=STRU%A(1)/2
        BR2(2,3)=STRU%A(3)/2
        BR2(3,1)=STRU%A(1)/2
        BR2(3,2)=STRU%A(2)/2
!       << conventional unit cell >>
        BR1(1,1)=STRU%A(1)
        BR1(2,2)=STRU%A(2)
        BR1(3,3)=STRU%A(3)
      ELSE IF(STRU%LATTIC(1:1).EQ.'B') THEN
! -------------------------------------------------------------------------
!       << body-centered lattice : B a b c * * *
!
!       a_1 = (-a/2, b/2, c/2)
!       a_2 = ( a/2,-b/2, c/2)
!       a_3 = ( a/2, b/2,-c/2)
!
!       orthorhombic, tetragonal, cubic
!
!       b_1 = (  0  , 1/b , 1/c )
!       b_2 = ( 1/a ,  0  , 1/c )
!       b_3 = ( 1/a , 1/b ,  0  )
! -------------------------------------------------------------------------
!       << primitive unit cell >>
        BR2(1,1)=-STRU%A(1)/2
        BR2(1,2)= STRU%A(2)/2
        BR2(1,3)= STRU%A(3)/2
        BR2(2,1)= STRU%A(1)/2
        BR2(2,2)=-STRU%A(2)/2
        BR2(2,3)= STRU%A(3)/2
        BR2(3,1)= STRU%A(1)/2
        BR2(3,2)= STRU%A(2)/2
        BR2(3,3)=-STRU%A(3)/2
!       << conventional unit cell >>
        BR1(1,1)=STRU%A(1)
        BR1(2,2)=STRU%A(2)
        BR1(3,3)=STRU%A(3)
     else if(stru%lattic(1:3)=='CXY' .and. stru%alpha(3)==90) then
! -------------------------------------------------------------------------
!       << base-centered (in the xy-plane) : CXY a b c alp bet 90 >>
!
!       Note: either alp or bet must be 90 degree
!
!       a_1 = (a*sin(bet)/2,-b*sin(alp)/2,a*cos(bet)/2-b*cos(alp)/2)
!       a_2 = (a*sin(bet)/2, b*sin(alp)/2,a*cos(bet)/2+b*cos(alp)/2)
!       a_3 = (    0       ,    0       ,            c            )
!
!       monoclinic, orthorhombic
!
!       b_1 = ( 1/(a*sin(bet)),-1/(b*sin(alp)), 0 )
!       b_2 = ( 1/(a*sin(bet)), 1/(b*sin(alp)), 0 )
!       b_3 = (-1/(c*tan(bet)),-1/(c*tan(alp)),1/c)
! -------------------------------------------------------------------------
        if(stru%alpha(1)/=90 .and. stru%alpha(2)/=90) &
             call croak('LATTIC NOT DEFINED: '//stru%lattic)
!       << primitive unit cell >>
        BR2(1,1)= STRU%A(1)/2*SIN(BETA)
        BR2(1,2)=-STRU%A(2)/2*SIN(ALPHA)
        BR2(1,3)= STRU%A(1)/2*COS(BETA) &
                - STRU%A(2)/2*COS(BETA)
        BR2(2,1)= STRU%A(1)/2*SIN(BETA)
        BR2(2,2)= STRU%A(2)/2*SIN(ALPHA)
        BR2(2,3)= STRU%A(1)/2*COS(BETA) &
                + STRU%A(2)/2*COS(ALPHA)
        BR2(3,3)= STRU%A(3)
!       << conventional unit cell >>
        BR1(1,1)=STRU%A(1)*SIN(BETA)
        BR1(1,3)=STRU%A(1)*COS(BETA)
        BR1(2,2)=STRU%A(2)*SIN(ALPHA)
        BR1(2,3)=STRU%A(2)*COS(ALPHA)
        BR1(3,3)=STRU%A(3)
     else if(stru%lattic(1:3)=='CYZ' .and. stru%alpha(1)==90) then
! -------------------------------------------------------------------------
!       << base-centered (in the yz-plane) : CYZ a b c 90 bet gam >>
!
!       Note: either bet or gam must be 90 degree
!
!       a_1 = (a*sin(bet)*sin(gam),a*sin(bet)*cos(gam),a*cos(bet))
!       a_2 = (        0        ,           b/2       ,  -c/2    )
!       a_3 = (        0        ,           b/2       ,   c/2    )
!
!       monoclinic, orthorhombic
!
!       b_1 = ( 1/(a*sin(bet)*sin(gam))               , 0 ,  0 )
!       b_2 = (-1/(b*tan(gam))+1/(c*tan(bet)*sin(gam)),1/b,-1/c)
!       b_3 = (-1/(b*tan(gam))-1/(c*tan(bet)*sin(gam)),1/b, 1/c)
! -------------------------------------------------------------------------
        if(stru%alpha(2)/=90 .and. stru%alpha(3)/=90) &
             call croak('LATTIC NOT DEFINED: '//stru%lattic)
!       << primitive unit cell >>
        BR2(1,1)= STRU%A(1)*SIN(BETA)*SIN(GAMMA)
        BR2(1,2)= STRU%A(1)*SIN(BETA)*COS(GAMMA)
        BR2(1,3)= STRU%A(1)*COS(BETA)
        BR2(2,2)= STRU%A(2)/2
        BR2(2,3)=-STRU%A(3)/2
        BR2(3,2)= STRU%A(2)/2
        BR2(3,3)= STRU%A(3)/2
!       << conventional unit cell >>
        BR1(1,1)=STRU%A(1)*SIN(BETA)*SIN(GAMMA)
        BR1(1,2)=STRU%A(1)*SIN(BETA)*COS(GAMMA)
        BR1(1,3)=STRU%A(1)*COS(BETA)
        BR1(2,2)=STRU%A(2)
        BR1(3,3)=STRU%A(3)
     else if(stru%lattic(1:3)=='CXZ' .and. stru%alpha(2)==90) then
! -------------------------------------------------------------------------
!       << base-centered (in the xz-plane) : CXZ a b c alp 90 gam >>
!
!       Note: either alp or gam must be 90 degree
!
!       a_1 = (a*sin(gam)/2,a*cos(gam)/2,  -c/2    )
!       a_2 = (    0       ,b*sin(alp)  ,b*cos(alp))
!       a_3 = (a*sin(gam)/2,a*cos(gam)/2,   c/2    )
!
!       monoclinic, orthorhombic
!
!       b_1 = ( 1/(a*sin(gam))         , 1/(c*tan(alp)),-1/c)
!       b_2 = (-1/(b*tan(gam)*sin(alp)), 1/(b*sin(alp)),  0 )
!       b_3 = ( 1/(a*sin(gam))         ,-1/(c*tan(alp)), 1/c)
! -------------------------------------------------------------------------
        if(stru%alpha(1)/=90 .and. stru%alpha(3)/=90) &
             call croak('LATTIC NOT DEFINED: '//stru%lattic)
!       << primitive unit cell >>
        BR2(1,1)= STRU%A(1)/2*SIN(GAMMA)
        BR2(1,2)= STRU%A(1)/2*COS(GAMMA)
        BR2(1,3)=-STRU%A(3)/2
        BR2(2,2)= STRU%A(2)*SIN(ALPHA)
        BR2(2,3)= STRU%A(2)*COS(ALPHA)
        BR2(3,1)= STRU%A(1)/2*SIN(GAMMA)
        BR2(3,2)= STRU%A(1)/2*COS(GAMMA)
        BR2(3,3)= STRU%A(3)/2
!       << conventional unit cell >>
        BR1(1,1)=STRU%A(1)*SIN(GAMMA)
        BR1(1,2)=STRU%A(1)*COS(GAMMA)
        BR1(2,2)=STRU%A(2)*SIN(ALPHA)
        BR1(2,3)=STRU%A(2)*COS(ALPHA)
        BR1(3,3)=STRU%A(3)
      ELSE
        STOP 'LATTIC NOT DEFINED'
      END IF
!
!     << find reciprocal lattice vectors (without a factor of 2pi) >>
      CALL GBASS(BR1,BR3,.FALSE.)
      CALL GBASS(BR2,BR4,.FALSE.)
!
      WRITE(unit_out,1000) (I,(BR1(I,J),J=1,3),I=1,3)
      WRITE(unit_out,1010) (I,(BR2(I,J),J=1,3),I=1,3)
      WRITE(unit_out,1020) (I,(BR3(I,J),J=1,3),I=1,3)
      WRITE(unit_out,1030) (I,(BR4(I,J),J=1,3),I=1,3)
 1000 FORMAT(/' REAL SPACE LATTICE VECTORS a1, a2, a3 (in bohr)' &
             /' -----------------------------------------------' &
             /' CONVENTIONAL UNIT CELL :'/(' a',I1,'     = ',3F12.7))
 1010 FORMAT(/' PRIMITIVE UNIT CELL :'/(' a',I1,'     = ',3F12.7))
 1020 FORMAT(/' RECIPROCAL LATTIC VECTORS b1, b2, b3 (in 1/bohr)' &
             /' ------------------------------------------------' &
             /' CONVENTIONAL UNIT CELL :'/(' b',I1,'/2pi = ',3F12.7))
 1030 FORMAT(/' PRIMITIVE UNIT CELL :'/(' b',I1,'/2pi = ',3F12.7))
      END


!!/---
!! Local Variables:
!! mode: f90
!! End:
!!\---
!!
!! Time-stamp: <2016-07-06 15:58:41 assman@faepop71.tu-graz.ac.at>
