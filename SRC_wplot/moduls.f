!!! wien2wannier/SRC_wplot/moduls.f
!!!
!!! $Id: moduls.f 166 2014-02-03 09:39:24Z assmann $

module param
  use const, only: BUFSZ

  implicit none

  private :: BUFSZ
  public

!!!/=== The following is from the old ‘param.inc_{r,c}’ ===================
!     >> Standard WIEN97 Parameters <<
!
!     Constant parameter definition
!
  INTEGER   NRAD, NSYM, NRF, NLOAT
  INTEGER   LOMAX, LMAX7
  PARAMETER          (NRAD=    981)
  PARAMETER          (NSYM=     48)
  PARAMETER          (LMAX7=     8)
  PARAMETER          (LOMAX=     3)
  PARAMETER          (NRF = 4)
  PARAMETER          (NLOAT =3)

!:17[
!     There are a couple of optional features in how to set-up the wave
!     function's augmentation coefficients. None of them are essential
!     for the wave functions themselves. However, to stay consistent with
!     LAPW1 and LAPW2 they have been added to LAPW7 as well. To allow
!     access to these optional features four logical control parameters
!     have been introduced which have to be set to .TRUE. to activate
!     the LAPW1/LAPW2 features. Any other choice of these parameters 
!     is possible as well, but the code has not been optimized for any
!     of these alternative choices. For more details see the "trick 17"
!     sections C:17[ ... C:17] in the code (Uwe Birkenheuer)
!      
  LOGICAL         KCONJG,USEROT,ADDLOC,MVATOM
  PARAMETER       (KCONJG=.TRUE.)
  PARAMETER       (USEROT=.TRUE.)
  PARAMETER       (ADDLOC=.TRUE.)
  PARAMETER       (MVATOM=.TRUE.)
!:17]
!!!!!\=== END param.inc ===================================================

  integer, parameter :: unit_in=5, unit_out=6, unit_grid=7, unit_struct=8
  integer, parameter :: unit_vector=10, unit_vsp=18, unit_inwf=31, unit_chk=32
  integer, parameter :: unit_psink=21, unit_psiarg=22, unit_def=1
  integer, parameter :: unit_tmp=90, unit_rot=5894

!!! The following is set by ‘wplot.f’ to be passed to ‘main.F’
  character(BUFSZ) :: vecfn, psinkfn, psiargfn, outfn, tmpfn
  integer          :: idx_wann=0, iproc=0
end module param

module struct
  use const, only: R8
  real(r8), allocatable :: POS(:,:),ZNUC(:),RMT(:)
  integer,allocatable :: MULT(:),IATNR(:)
end module struct

module radgrd 
  use const, only: R8
  real(r8),allocatable :: RM(:,:),RNOT(:),DX(:)
  integer,allocatable :: JRI(:)
end module radgrd

module lolog 
  integer  NLO
  logical,allocatable :: LAPW(:,:)
  integer, allocatable :: ILO(:,:)
end module lolog

module loabc
  use const, only: R8
  real(r8),allocatable :: ALO(:,:,:,:)  
end module loabc

module atspdt
  use const, only: R8
  real(r8),allocatable :: P(:,:,:),DP(:,:,:) 
end module atspdt

module radfu
  use const, only: R8
  real(r8),allocatable ::  RRAD(:,:,:,:)
end module radfu

module bessfu
  use const, only: R8
  real(r8),allocatable ::  FJ(:,:,:),DFJ(:,:,:),RAD(:)
  integer, allocatable :: IRAD(:)
end module bessfu

module work
  use const, only: C16
  complex(c16),allocatable :: aug(:,:,:)
end module work

module grid
  use const, only: R8
  real(r8),allocatable :: rgrid(:,:)
  integer,allocatable :: ireg(:),ilat(:,:),iri(:)
  integer npg
end module grid

MODULE ams
  use const, only: R8
  REAL(R8)   :: atom_mass(103)

CONTAINS
  SUBROUTINE init_ams
    REAL(R8) :: atom_mass(103)
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

module latt
  use const, only: R8

  implicit none

  private R8
  public

  real(R8) :: VUC, BR1(3,3), BR2(3,3), BR3(3,3), BR4(3,3)
end module latt

module sym2
  use const, only: R8
  use param, only: NSYM

  implicit none

  private :: R8, NSYM
  public

  real(R8) :: TAU(3,NSYM)
  integer  :: IMAT(3,3,NSYM), IORD
end module sym2

!!/---
!! Local Variables:
!! mode: f90
!! End:
!!\---
