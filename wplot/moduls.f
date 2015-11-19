      module struct
      real*8, allocatable :: POS(:,:),ZNUC(:),RMT(:)
      integer,allocatable :: MULT(:),IATNR(:)
      end module struct
!
      module radgrd 
      real*8,allocatable :: RM(:,:),RNOT(:),DX(:)
      integer,allocatable :: JRI(:)
      end module radgrd
!
      module lolog 
      integer  NLO
      logical,allocatable :: LAPW(:,:)
      integer, allocatable :: ILO(:,:)
      end module lolog
!
      module loabc
      real*8,allocatable :: ALO(:,:,:,:)  
      end module loabc
!
      module atspdt
      real*8,allocatable :: P(:,:,:),DP(:,:,:) 
      end module atspdt
!
      module radfu
      real*8,allocatable ::  RRAD(:,:,:,:)
      end module radfu
!
      module bessfu
      real*8,allocatable ::  FJ(:,:,:),DFJ(:,:,:),RAD(:)
      integer, allocatable :: IRAD(:)
      end module bessfu
!
      module work
      complex*16,allocatable :: aug(:,:,:)
      end module work
!
      module grid
      real*8,allocatable :: rgrid(:,:)
      integer,allocatable :: ireg(:),ilat(:,:),iri(:)
      integer npg
      end module grid

MODULE ams
  REAL*8          :: atom_mass(103)
    
 CONTAINS
  SUBROUTINE init_ams
  REAL*8          :: atom_mass(103)
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

!
