      SUBROUTINE LATGEN(LATTIC,A,ALP,BET,GAM,ORTHO,PRIM)
      IMPLICIT    REAL*8 (A-H,O-Z)
      CHARACTER*4 LATTIC
      DIMENSION   A(3)
      LOGICAL     ORTHO, PRIM
      COMMON /LATT/ VUC,BR1(3,3),BR2(3,3),BR3(3,3),BR4(3,3)
!
! << Input >>
! LATTIC      -- the type of the Bravais lattice
! A(1:3)      -- the lattice constants
! ALP,BET,GAM -- the unit cell angles (in degree)
!
! << Output >>
! ORTHO       -- .TRUE. if the lattice is an orthogonal one
! PRIM        -- .TRUE. if the lattice is a primitive one
!
! << Output (to COMMON /LATT/) >>
! VUC      -- the volume of the primitive unit cell
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
      ALPHA = ALP * 3.14159265358979324D0 / 180.0D0
      BETA  = BET * 3.14159265358979324D0 / 180.0D0
      GAMMA = GAM * 3.14159265358979324D0 / 180.0D0
      DO 10 I=1,3
      DO 10 J=1,3
      BR2(I,J)=0.0D0
   10 BR1(I,J)=0.0D0
!
      IF(LATTIC(1:1).EQ.'P')THEN
! -------------------------------------------------------------------------
!       << primitive lattice : P a b c alp bet gam >>
!
!       a_1 = a * (sin(bet)*sin(phi),sin(bet)*cos(phi),cos(bet))
!       a_2 = b * (        0        ,      sin(alp)   ,cos(alp))
!       a_3 = c * (        0        ,        0        ,  1     )
!       with
!       cos(phi) := ( cos(gam) - cos(alp)*cos(bet) ) / sin(alp)*sin(bet)
!
!       VUC = a*b*c * sin(alp)*sin(bet)*sin(phi)
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
        BR2(1,1)=A(1)*SIN(BETA)*SIN(PHI)
        BR2(1,2)=A(1)*SIN(BETA)*COS(PHI)
        BR2(1,3)=A(1)*COS(BETA)
        BR2(2,2)=A(2)*SIN(ALPHA)
        BR2(2,3)=A(2)*COS(ALPHA)
        BR2(3,3)=A(3)
!       << conventional unit cell >>
        DO 11 I=1,3
        DO 11 J=1,3
   11   BR1(I,J) = BR2(I,J)
!       << further settings >>
        VUC   = SIN(ALPHA)*SIN(BETA)*SIN(PHI) * A(1)*A(2)*A(3)
        PRIM  = .TRUE.
        ORTHO = ALP.EQ.90.0D0 .AND. BET.EQ.90.0D0 .AND. GAM.EQ.90.0D0
      ELSE IF(LATTIC(1:1).EQ.'H') THEN
! -------------------------------------------------------------------------
!       << hexagonal lattice : H a * c * * * >>
!
!       a_1 = a * (sqrt(3)/2,-1/2,0)
!       a_2 = a * (     0   ,  1 ,0)
!       a_3 = c * (     0   ,  0 ,1)
!
!       VUC = a^2*c * sqrt(3)/2
!
!       this setting corresponds to the P a a c 90 90 120 setting
!
!       b_1 = ( 1/sqrt(3),  -1 ,   1 ) / a
!       b_2 = ( 1/sqrt(3),   1 ,   1 ) / a
!       b_3 = (-2/sqrt(3),   0 ,   1 ) / c
! -------------------------------------------------------------------------
!       << primitive unit cell >>
        BR2(1,1)= A(1)*SQRT(3.0D0)/2.0D0
        BR2(1,2)=-A(1)*0.5D0
        BR2(2,2)= A(1)
        BR2(3,3)= A(3)
!       << conventional unit cell >>
        DO 12 I=1,3
        DO 12 J=1,3
   12   BR1(I,J) = BR2(I,J)
!       << further settings >>
        VUC   = A(1)**2*A(3) * SQRT(3.0D0)/2.0D0
        PRIM  = .TRUE.
        ORTHO = .FALSE.
      ELSE IF(LATTIC(1:1).EQ.'T'.OR.LATTIC(1:1).EQ.'R') THEN
! -------------------------------------------------------------------------
!       << rhombohedral or trigonal lattice : R a * c * * * >>
!
!       a_1 = ( a/(2*sqrt(3)),-a/2,c/3)
!       a_2 = ( a/(2*sqrt(3)), a/2,c/3)
!       a_3 = (-a/   sqrt(3) ,  0 ,c/3)
!
!       VUC = a^2*c / (2*sqrt(3))
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
        BR2(1,1)= A(1)*0.5D0/SQRT(3.0D0)
        BR2(1,2)=-A(1)*0.5D0
        BR2(1,3)= A(3)/3.0D0
        BR2(2,1)= A(1)*0.5D0/SQRT(3.0D0)
        BR2(2,2)= A(1)*0.5D0
        BR2(2,3)= A(3)/3.0D0
        BR2(3,1)=-A(1)/SQRT(3.0D0)
        BR2(3,3)= A(3)/3.0D0
!       << conventional unit cell >>
        DO 13 I=1,3
        DO 13 J=1,3
   13   BR1(I,J) = BR2(I,J)
!       << further settings >>
        VUC   = A(1)**2*A(3) * 0.5D0/SQRT(3.0D0)
        PRIM  = .TRUE.
        ORTHO = .FALSE.
      ELSE IF(LATTIC(1:1).EQ.'F') THEN
! -------------------------------------------------------------------------
!       << face-centered lattice : F a b c * * * 
!
!       a_1 = ( 0 ,b/2,c/2)
!       a_2 = (a/2, 0 ,c/2)
!       a_3 = (a/2,b/2, 0 )
!
!       VUC = a*b*c / 4
!
!       orthorhombic, cubic
!
!       b_1 = ( -1/a ,  1/b ,  1/c )
!       b_2 = (  1/a , -1/b ,  1/c )
!       b_3 = (  1/a ,  1/b , -1/c )
! -------------------------------------------------------------------------
!       << primitive unit cell >>
        BR2(1,2)=A(2)*0.5D0
        BR2(1,3)=A(3)*0.5D0
        BR2(2,1)=A(1)*0.5D0
        BR2(2,3)=A(3)*0.5D0
        BR2(3,1)=A(1)*0.5D0
        BR2(3,2)=A(2)*0.5D0
!       << conventional unit cell >>
        BR1(1,1)=A(1)
        BR1(2,2)=A(2)
        BR1(3,3)=A(3)
!       << further settings >>
        VUC   = A(1)*A(2)*A(3) / 4.0D0
        PRIM  = .FALSE.
        ORTHO = .TRUE.
      ELSE IF(LATTIC(1:1).EQ.'B') THEN
! -------------------------------------------------------------------------
!       << body-centered lattice : B a b c * * *
!
!       a_1 = (-a/2, b/2, c/2)
!       a_2 = ( a/2,-b/2, c/2)
!       a_3 = ( a/2, b/2,-c/2)
!
!       VUC = a*b*c / 2
!
!       orthorhombic, tetragonal, cubic
!
!       b_1 = (  0  , 1/b , 1/c )
!       b_2 = ( 1/a ,  0  , 1/c )
!       b_3 = ( 1/a , 1/b ,  0  )
! -------------------------------------------------------------------------
!       << primitive unit cell >>
        BR2(1,1)=-A(1)*0.5D0
        BR2(1,2)= A(2)*0.5D0
        BR2(1,3)= A(3)*0.5D0
        BR2(2,1)= A(1)*0.5D0
        BR2(2,2)=-A(2)*0.5D0
        BR2(2,3)= A(3)*0.5D0
        BR2(3,1)= A(1)*0.5D0
        BR2(3,2)= A(2)*0.5D0
        BR2(3,3)=-A(3)*0.5D0
!       << conventional unit cell >>
        BR1(1,1)=A(1)
        BR1(2,2)=A(2)
        BR1(3,3)=A(3)
!       << further settings >>
        VUC   = A(1)*A(2)*A(3) / 2.0D0
        PRIM  = .FALSE.
        ORTHO = .TRUE.
      ELSE IF(LATTIC(1:3).EQ.'CXY' .AND. GAM.EQ.90.0D0) THEN
! -------------------------------------------------------------------------
!       << base-centered (in the xy-plane) : CXY a b c alp bet 90 >>
!
!       Note: either alp or bet must be 90 degree
!
!       a_1 = (a*sin(bet)/2,-b*sin(alp)/2,a*cos(bet)/2-b*cos(alp)/2)
!       a_2 = (a*sin(bet)/2, b*sin(alp)/2,a*cos(bet)/2+b*cos(alp)/2)
!       a_3 = (    0       ,    0       ,            c            )
!
!       VUC = a*b*c * sin(alp)*sin(bet)/2
!       
!       monoclinic, orthorhombic
!       
!       b_1 = ( 1/(a*sin(bet)),-1/(b*sin(alp)), 0 )
!       b_2 = ( 1/(a*sin(bet)), 1/(b*sin(alp)), 0 )
!       b_3 = (-1/(c*tan(bet)),-1/(c*tan(alp)),1/c)
! -------------------------------------------------------------------------
        IF(ALP.NE.90.0D0 .AND. BET.NE.90.0D0) STOP 'LATTIC NOT DEFINED'
!       << primitive unit cell >>
        BR2(1,1)= A(1)*0.5D0*SIN(BETA)
        BR2(1,2)=-A(2)*0.5D0*SIN(ALPHA)
        BR2(1,3)= A(1)*0.5D0*COS(BETA) &
                - A(2)*0.5D0*COS(BETA)
        BR2(2,1)= A(1)*0.5D0*SIN(BETA)
        BR2(2,2)= A(2)*0.5D0*SIN(ALPHA)
        BR2(2,3)= A(1)*0.5D0*COS(BETA) &
                + A(2)*0.5D0*COS(ALPHA)
        BR2(3,3)= A(3)
!       << conventional unit cell >>
        BR1(1,1)=A(1)*SIN(BETA)
        BR1(1,3)=A(1)*COS(BETA)
        BR1(2,2)=A(2)*SIN(ALPHA)
        BR1(2,3)=A(2)*COS(ALPHA)
        BR1(3,3)=A(3)
!       << further settings >>
        VUC   = A(1)*A(2)*A(3) * SIN(ALPHA)*SIN(BETA)/2.0D0
        PRIM  = .FALSE.
        ORTHO = ALP.EQ.90.0D0 .AND. BET.EQ.90.0D0
      ELSE IF(LATTIC(1:3).EQ.'CYZ' .AND. ALP.EQ.90.0D0) THEN
! -------------------------------------------------------------------------
!       << base-centered (in the yz-plane) : CYZ a b c 90 bet gam >>
!
!       Note: either bet or gam must be 90 degree
!
!       a_1 = (a*sin(bet)*sin(gam),a*sin(bet)*cos(gam),a*cos(bet))
!       a_2 = (        0        ,           b/2       ,  -c/2    )
!       a_3 = (        0        ,           b/2       ,   c/2    )
!
!       VUC = a*b*c * sin(bet)*sin(gam)/2
!       
!       monoclinic, orthorhombic
!       
!       b_1 = ( 1/(a*sin(bet)*sin(gam))               , 0 ,  0 )
!       b_2 = (-1/(b*tan(gam))+1/(c*tan(bet)*sin(gam)),1/b,-1/c)
!       b_3 = (-1/(b*tan(gam))-1/(c*tan(bet)*sin(gam)),1/b, 1/c)
! -------------------------------------------------------------------------
        IF(BET.NE.90.0D0 .AND. GAM.NE.90.0D0) STOP 'LATTIC NOT DEFINED'
!       << primitive unit cell >>
        BR2(1,1)= A(1)*SIN(BETA)*SIN(GAMMA)
        BR2(1,2)= A(1)*SIN(BETA)*COS(GAMMA)
        BR2(1,3)= A(1)*COS(BETA)
        BR2(2,2)= A(2)*0.5D0
        BR2(2,3)=-A(3)*0.5D0
        BR2(3,2)= A(2)*0.5D0
        BR2(3,3)= A(3)*0.5D0
!       << conventional unit cell >>
        BR1(1,1)=A(1)*SIN(BETA)*SIN(GAMMA)
        BR1(1,2)=A(1)*SIN(BETA)*COS(GAMMA)
        BR1(1,3)=A(1)*COS(BETA)
        BR1(2,2)=A(2)
        BR1(3,3)=A(3)
!       << further settings >>
        VUC   = A(1)*A(2)*A(3) * SIN(BETA)*SIN(GAMMA)/2.0D0
        PRIM  = .FALSE.
        ORTHO = BET.EQ.90.0D0 .AND. GAM.EQ.90.0D0
      ELSE IF(LATTIC(1:3).EQ.'CXZ' .AND. BET.EQ.90.0D0) THEN
! -------------------------------------------------------------------------
!       << base-centered (in the xz-plane) : CXZ a b c alp 90 gam >>
!
!       Note: either alp or gam must be 90 degree
!
!       a_1 = (a*sin(gam)/2,a*cos(gam)/2,  -c/2    )
!       a_2 = (    0       ,b*sin(alp)  ,b*cos(alp))
!       a_3 = (a*sin(gam)/2,a*cos(gam)/2,   c/2    )
!
!       VUC = a*b*c * sin(alp)*sin(gam)/2
!       
!       monoclinic, orthorhombic
!
!       b_1 = ( 1/(a*sin(gam))         , 1/(c*tan(alp)),-1/c)
!       b_2 = (-1/(b*tan(gam)*sin(alp)), 1/(b*sin(alp)),  0 )
!       b_3 = ( 1/(a*sin(gam))         ,-1/(c*tan(alp)), 1/c)
! -------------------------------------------------------------------------
        IF(ALP.NE.90.0D0 .AND. GAM.NE.90.0D0) STOP 'LATTIC NOT DEFINED'
!       << primitive unit cell >>
        BR2(1,1)= A(1)*0.5D0*SIN(GAMMA)
        BR2(1,2)= A(1)*0.5D0*COS(GAMMA)
        BR2(1,3)=-A(3)*0.5D0
        BR2(2,2)= A(2)*SIN(ALPHA)
        BR2(2,3)= A(2)*COS(ALPHA)
        BR2(3,1)= A(1)*0.5D0*SIN(GAMMA)
        BR2(3,2)= A(1)*0.5D0*COS(GAMMA)
        BR2(3,3)= A(3)*0.5D0
!       << conventional unit cell >>
        BR1(1,1)=A(1)*SIN(GAMMA)
        BR1(1,2)=A(1)*COS(GAMMA)
        BR1(2,2)=A(2)*SIN(ALPHA)
        BR1(2,3)=A(2)*COS(ALPHA)
        BR1(3,3)=A(3)
!       << further settings >>
        VUC   = A(1)*A(2)*A(3) * SIN(ALPHA)*SIN(GAMMA)/2.0D0
        PRIM  = .FALSE.
        ORTHO = ALP.EQ.90.0D0 .AND. GAM.EQ.90.0D0
      ELSE
        STOP 'LATTIC NOT DEFINED'
      END IF
!
!     << find reciprocal lattice vectors (without a factor of 2pi) >>
      CALL GBASS(BR1,BR3,.FALSE.)
      CALL GBASS(BR2,BR4,.FALSE.)
!
      WRITE(6,1000) (I,(BR1(I,J),J=1,3),I=1,3)
      WRITE(6,1010) (I,(BR2(I,J),J=1,3),I=1,3)
      WRITE(6,1020) (I,(BR3(I,J),J=1,3),I=1,3)
      WRITE(6,1030) (I,(BR4(I,J),J=1,3),I=1,3)
 1000 FORMAT(/' REAL SPACE LATTICE VECTORS a1, a2, a3 (in bohr)' &
             /' -----------------------------------------------' &
             /' CONVENTIONAL UNIT CELL :'/(' a',I1,'     = ',3F12.7))
 1010 FORMAT(/' PRIMITIVE UNIT CELL :'/(' a',I1,'     = ',3F12.7))
 1020 FORMAT(/' RECIPROCAL LATTIC VECTORS b1, b2, b3 (in 1/bohr)' &
             /' ------------------------------------------------' &
             /' CONVENTIONAL UNIT CELL :'/(' b',I1,'/2pi = ',3F12.7))
 1030 FORMAT(/' PRIMITIVE UNIT CELL :'/(' b',I1,'/2pi = ',3F12.7))
      END
