!!! wien2wannier/SRC_w2w/modules.F

MODULE param
  use const, only: R8, C16
  private :: R8, C16

  public

  character(*), parameter, private :: rev_str="$version: v1.0.0-120-g28e9feb$"
  character(*), parameter, public  :: &
       wien2wannier_version = rev_str(11 : len (rev_str)-1)

  integer, parameter :: unit_in=5, unit_out=6, unit_amn=7, unit_mmn=8
  integer, parameter :: unit_vector=9, unit_nnkp=10, unit_eig=12, unit_vsp=18
  integer, parameter :: unit_struct=20, unit_ene=50, unit_def=1, unit_fermi=51

  integer            LMAX2, LOMAX
  integer            NCOM, NGAU
  integer            NRAD!, NPAIR, NBANDS
  integer            IBLOCK
  integer            LXDOS, NLOAT, NRF
  real(r8)             clight
  PARAMETER (IBLOCK= 128)
  !.....Optimize IBLOCK and IBLCK for your hardware (32-255)
  PARAMETER (LMAX2=  5)
  PARAMETER (LOMAX=   3)
  INTEGER        :: NATO =    0
  PARAMETER (NCOM=   49)
  ! for ncom parameter check format 1003 in l2main.frc
  INTEGER        :: NDIF =    0
  PARAMETER (NGAU= 1500)
  INTEGER        :: NMAT =    0
  PARAMETER (NRAD=  881)
  INTEGER        :: NSYM =    0
  INTEGER        :: NUME =    0
  ! for x-dos set lxdos to 3
  parameter (lxdos= 3)
  parameter (nloat= 3)
  parameter (nrf=4)
  PARAMETER (CLIGHT=137.0359895D0)
END MODULE param

module assleg
  !      DOUBLE PRECISION   YR(N,MAXDIM)
  !      COMMON  /ASSLEG/   YR
  !      SAVE    /ASSLEG/
  DOUBLE PRECISION, allocatable :: YR(:,:)
  INTEGER N, MAXDIM
  PARAMETER (N=6, MAXDIM=81)
contains
  SUBROUTINE init_assleg
    IMPLICIT NONE
    !          integer N,MAXDIM
    allocate( YR(N, MAXDIM) )
    YR=0.0D0
  END SUBROUTINE init_assleg
end module assleg


MODULE struct
  USE param

  use const, only: R8
  private :: R8

  public

  LOGICAL                    :: rel
  REAL(R8)                   :: AA,BB,CC,VOL,pia(3),alpha(3)
  REAL(R8),ALLOCATABLE       :: R0(:),DX(:),RMT(:),zz(:),rotloc(:,:,:),v(:)
  REAL(R8),ALLOCATABLE       :: tau(:,:)
  REAL(R8),ALLOCATABLE       :: pos(:,:)
  REAL(R8),ALLOCATABLE       :: rotij(:,:,:),tauij(:,:)
  CHARACTER(4)               :: lattic,irel,cform
  CHARACTER(80)              :: title
  CHARACTER(10), ALLOCATABLE :: aname(:)
  INTEGER                    :: nat,iord
  INTEGER,ALLOCATABLE        :: mult(:),jrj(:)
  INTEGER,ALLOCATABLE        :: iz(:,:,:)

 CONTAINS
  SUBROUTINE init_struct
    USE reallocate

    use const, only: R8

    IMPLICIT NONE

    INTEGER                :: ios,inum,isplit,iatnr
    REAL(R8)               :: test,ninety
    INTEGER                :: index,i,j,j1,j2,m,jatom

    test=1.D-5
    ninety=90.0D0

    read (unit_struct,1000) title
    read (unit_struct,1010)  lattic,nat,cform,irel
    REL=.TRUE.
    IF(IREL.EQ.'NREL') REL=.FALSE.
    ALLOCATE(aname(nat),mult(nat),jrj(nat),r0(nat),dx(nat),rmt(nat))
    allocate(zz(nat),rotloc(3,3,nat),v(nat))
    v=0.0d0
    ALLOCATE (pos(3,48*nat))
    read (unit_struct,1020) aa,bb,cc,alpha(1),alpha(2),alpha(3)
    IF(ABS(ALPHA(1)).LT.test) ALPHA(1)=ninety
    IF(ABS(ALPHA(2)).LT.test) ALPHA(2)=ninety
    IF(ABS(ALPHA(3)).LT.test) ALPHA(3)=ninety
    INDEX=0
    DO jatom=1,NAT
     INDEX=INDEX+1
     read(unit_struct,1030,iostat=ios) &
          iatnr, (pos(j,index),j=1,3), mult(jatom), isplit
     IF(ios /= 0) THEN
      WRITE(unit_out,*) iatnr,(pos(j,index),j=1,3),mult(jatom)
      WRITE(unit_out,*) 'ERROR IN STRUCT FILE READ'
      STOP
      ENDIF
      IF (mult(jatom) .EQ. 0) THEN
       WRITE (unit_out,6000) jatom, index, mult(jatom)
       STOP
       ENDIF
       DO m=1,mult(jatom)-1
        index=index+1
        READ(unit_struct,1031)iatnr,(pos(j,index),j=1,3)
       ENDDO
       READ(unit_struct,1050) aname(jatom),jrj(jatom),r0(jatom),rmt(jatom), &
            zz(jatom)
       dx(jatom)=LOG(rmt(jatom)/r0(jatom)) / (jrj(jatom)-1)
       rmt(jatom)=r0(jatom)*EXP( dx(jatom)*(jrj(jatom)-1) )
       READ(unit_struct,1051) ((rotloc(i,j,jatom),i=1,3),j=1,3)
    ENDDO
    ndif=index
    call realloc(pos, (/ 3, ndif /))
    ALLOCATE(rotij(3,3,ndif),tauij(3,ndif))
    READ(unit_struct,1151)iord
    ALLOCATE(iz(3,3,iord),tau(3,iord))
    DO j=1,iord
       READ(unit_struct,1101)((iz(j1,j2,j),j1=1,3),tau(j2,j),j2=1,3),inum
    ENDDO

1000 FORMAT(A80)
1010 FORMAT(A4,23X,I3,1x,a4,/,13X,A4,18X,A4)
1020 FORMAT(6F10.7,10X,F10.7)
1030 FORMAT(4X,I4,4X,F10.7,3X,F10.7,3X,F10.7,/,15X,I2,17X,I2)
1031 FORMAT(4X,I4,4X,F10.7,3X,F10.7,3X,F10.7)
1050 FORMAT(A10,5X,I5,5X,F10.9,5X,F10.5,5X,F10.5)
1051 FORMAT(20X,3F10.8)
1101 FORMAT(3(3I2,F11.8/),I8)
1151 FORMAT(I4)
6000 FORMAT(///,3X,'ERROR IN LAPW0 : MULT(JATOM)=0 ...', &
          /, 20X,'JATOM=',I3,3X,'INDEX=',I3,3X,'MULT=',I3)
  END SUBROUTINE init_struct

END MODULE struct

MODULE bessel
  use const
  REAL(R8), ALLOCATABLE  ::rj(:,:),ri_mat(:,:)

CONTAINS
  SUBROUTINE init_bessel(LMAX2,LJMAX,NRAD,NRF)
    L=0
    DO L1=0,LMAX2
       DO L2=0,LMAX2
          DO LJ=0,LJMAX
             IF (MOD((L1+L2+LJ),2) .EQ. 1) cycle
             IF ((L1+L2-LJ) .LT. 0)        cycle
             IF ((L1-L2+LJ) .LT. 0)        cycle
             IF ((-L1+L2+LJ) .LT. 0)       cycle
             L=L+1
          END DO
       END DO
    END DO
    ALLOCATE(rj(0:LJMAX+1,NRAD),ri_mat(NRF*NRF,L))
  END SUBROUTINE init_bessel
END MODULE bessel



MODULE xa
  use const

  COMPLEX(C16),ALLOCATABLE :: phs(:)
  REAL(R8),ALLOCATABLE     :: fj(:,:),dfj(:,:),r(:)
  REAL(R8)                 :: bk(3),bkrot(3),bkrloc(3)

 CONTAINS
  SUBROUTINE init_xa(LMAX2,NMAT,NRAD,NB)
    ALLOCATE(phs(nb),fj(0:lmax2,nmat),dfj(0:lmax2,nmat))
    ALLOCATE(r(nrad))
  END SUBROUTINE init_xa
END MODULE xa

MODULE xa3
  use const
#ifndef _COMPLEX_
  REAL(R8),ALLOCATABLE  :: a(:,:,:),apw1(:,:),apw2(:,:)
#else
  COMPLEX(C16),ALLOCATABLE  :: a(:,:,:),apw1(:,:),apw2(:,:)
#endif
  INTEGER,ALLOCATABLE  :: size(:),GX(:,:),GY(:,:),GZ(:,:)
  REAL(R8),ALLOCATABLE    :: BKX(:),BKY(:),BKZ(:)
  REAL(R8),ALLOCATABLE    :: XK(:),YK(:),ZK(:)

CONTAINS
  SUBROUTINE init_xa3(nb,nmat,num_kpts)
    ALLOCATE(apw1(Nb,Nmat), apw2(Nb,Nmat), a(Nmat,Nb,num_kpts))
    ALLOCATE(BKX(Nmat), BKY(Nmat), BKZ(Nmat))
    ALLOCATE(GX(Nmat,num_kpts), GY(Nmat,num_kpts), GZ(Nmat,num_kpts))
    ALLOCATE(XK(num_kpts), YK(num_kpts), ZK(num_kpts), size(num_kpts))
    GX = 0; GY = 0; GZ = 0
  END SUBROUTINE init_xa3
END MODULE xa3

MODULE amn_mmn
  USE param

  use const, only: C16
  private :: C16

  COMPLEX(C16),ALLOCATABLE     :: overlap(:,:,:),c(:,:,:)

CONTAINS
  SUBROUTINE init_mmn(nbands,npair,nproj)
    USE param
    INTEGER    :: nbands,npair,nproj
    ALLOCATE(overlap(nbands,nbands,npair),c(nproj,(LMAX2+1)*(LMAX2+1),ndif))
  END SUBROUTINE init_mmn
END MODULE amn_mmn

MODULE pairs
  INTEGER,ALLOCATABLE       :: KP(:),KPB(:),BQX(:),BQY(:),BQZ(:),BQX1(:),BQY1(:),BQZ1(:)

CONTAINS
  SUBROUTINE  init_pairs(npair)
    INTEGER  :: npair
    ALLOCATE(kp(npair),kpb(npair),bqx(npair),bqy(npair),bqz(npair),bqx1(npair),bqy1(npair),bqz1(npair))
  END SUBROUTINE init_pairs
END MODULE pairs


MODULE ams
  use const, only: R8

  REAL(R8)          :: atom_mass(103)

CONTAINS
  SUBROUTINE init_ams
    REAL(R8)          :: atom_mass(103)
    DATA atom_mass /1.,4.,6.9,9.,10.8,12.,14.,16.,19.,20.2, &
         23.,24.3,27.,28.1,31.,32.,35.4,40.,39.1,40.,45., &
         47.9,50.9,52.,54.9,55.8,58.9,58.7,63.5,65.4,69.7, &
         72.6,74.9,79.,79.9,83.8,85.5,87.6,88.9,91.2,92.9, &
         95.9,98.,101.1,102.9,106.4,107.9,112.4,114.8, &
         118.7,121.8,127.6,126.9,131.3,132.9,137.3,138.9,140.1, &
         140.9,144.2,145.,150.4,152.,157.3,158.9,162.5, &
         164.9,167.3,168.9,173.,175.,178.5,180.9,183.8,186.2, &
         190.2,192.2,195.1,197.,200.6,204.4,207.2,209.,209., &
         210.,222.,223.,226.,227.,232.,231.,238.,237.,244.,243., &
         247.,247.,251.,252.,257.,258.,259.,262./
  END SUBROUTINE init_ams
END MODULE ams

module lolog
  use param, only: lmax2, lomax

  implicit none

  integer :: nlo, nlov, nlon, n_rad(0:lmax2), ilo(0:lomax)
  logical :: loor(0:lomax), lapw(0:lmax2)
end module lolog

module gener
  use const, only: R8

  ! transformation between u.c. and cartesian coordinates
  real(R8) :: BR1(3,3), BR2(3,3)
end module gener

!!/---
!! Local Variables:
!! mode: f90
!! End:
!!\---
!!

!! Time-stamp: <2016-04-08 17:10:31 assman@faepop36.tu-graz.ac.at>
