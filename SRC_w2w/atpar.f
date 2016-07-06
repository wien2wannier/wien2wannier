!!! wien2wannier/SRC_w2w/atpar.f

subroutine ATPAR (JATOM, itape, jtape)
  ! calculate radial functions for atoms JATOM
  use param,  only: unit_out, clight, Nrad, Nloat, lomax, Lmax2
  use struct, only: JRJ, mult, Nat, aname, R0, dx, zz, rel
  use lolog,  only: nlo,nlov,nlon,loor,ilo,lapw,n_rad
  use atspdt, only: P, DP
  use const,  only: R8
  use uhelp,  only: A, B
  use radfu,  only: RF1, RF2

  implicit none

  integer, intent(in) :: jatom, itape, jtape

  real(R8) :: VR(Nrad), AE(Nrad), BE(Nrad)
  logical  :: rlo(1:nloat, 0:lomax)
  real(r8) :: emist(0:lomax,nloat),E(0:LMAX2),elo(0:LOMAX,nloat),pei(0:lmax2)
  integer  :: imax,irf, jlo, kappa, i,j,k,l, node,nodes,nodel, m
  real(R8) :: dele,delei, fl, ei,e1, uvb,duvb,uv,duv,uve,duve, ovlp, trx
  real(R8) :: try, r_m, pi12lo, pe12lo, cross

2022 FORMAT(3X,4E19.12)

  !.....READ TOTAL SPHERICAL POTENTIAL V(0,0) OF TAPEjtape=VSP
  !     NORM OF V(0,0)=V00(R)*R/SQRT(4.D0*PI)

  READ(jtape,1980)
  READ(jtape,2000)
  READ(jtape,2031)
  READ(jtape,2022)(VR(J),J=1,JRJ(JATOM))
  READ(jtape,2031)
  READ(jtape,2030)

  do J=1,JRJ(JATOM)
     VR(J)=VR(J)/2.0D0
  end do

  write(unit_out,*)'ATPAR'
  nlo=0
  nlov=0
  nlon=0
  ilo=0
  n_rad=2

  atoms: do I=1,JATOM
     READ(itape)E
     READ(itape)elo
     IF(i.EQ.jatom) THEN
        DO l=0,lmax2
           lapw(l)=.TRUE.
           IF(e(l).GT.150.) THEN
              e(l)=e(l)-200.d+0
              lapw(l)=.FALSE.
           ENDIF
        ENDDO
     ENDIF
     LOs: do l = 0,lomax
        loor(l)=.FALSE.
        do k=1,nloat
           rlo(k,l)=.FALSE.
           IF (i.EQ.jatom) THEN
              IF (elo(l,k).LT.(995.d+0)) THEN
                 ilo(l)=ilo(l)+1
                 nlo=nlo+((2*l+1))*mult(i)
                 IF(.NOT.lapw(l).AND.k.EQ.1) cycle
                 IF(k.EQ.nloat) THEN
                    rlo(ilo(l),l)=.TRUE.
                    cycle
                 ENDIF
                 loor(l)=.TRUE.
              ENDIF
           ELSE
              IF (elo(l,k).LT.(995.d+0)) nlov=nlov+((2*l+1))*mult(i)
           ENDIF
        end do
     end do LOs
  end do atoms

  if(jatom /= nat) then
     do I=JATOM+1,NAT
        READ(itape) EMIST
        READ(itape) EMIST
        do l=0,lomax
           do k=1,nloat
              if (emist(l,k).lt.(995.0D+0))  &
                   nlon=nlon+((2*l+1))*mult(i)
           end do
        end do
     end do
  end if

  WRITE(unit_out,7) ANAME(JATOM)
  WRITE(unit_out,5) E
  WRITE(unit_out,14)

  lloop: DO l=0,LMAX2
     DELE=2.0D-3
     DELEI=0.25D0/DELE
     FL=L
     EI=E(l)/2.0d0
     !     CALCULATE ENERGY-DERIVATIVE BY FINITE DIFFERENCE
     !     DELE IS THE UPWARD AND DOWNWARD ENERGY SHIFT IN HARTREES

     E1=EI-DELE
     CALL OUTWIN(REL,VR,R0(JATOM),DX(JATOM),JRJ(JATOM),E1,            &
          FL,UVB,DUVB,NODEL,ZZ(jatom))
     CALL RINT13(A,B,A,B,OVLP,JATOM)
     TRX=1.0D0/SQRT(OVLP)
     IMAX=JRJ(JATOM)
     DO M=1,IMAX
        AE(M)=TRX*A(M)
        BE(M)=TRX*B(M)
     end do
     UVB=TRX*UVB
     DUVB=TRX*DUVB
     E1=EI+DELE
     CALL OUTWIN(REL,VR,R0(JATOM),DX(JATOM),JRJ(JATOM),E1,            &
          FL,UVE,DUVE,NODE,ZZ(jatom))
     CALL RINT13(A,B,A,B,OVLP,JATOM)
     TRX=1.0d0/SQRT(OVLP)
     UVE=DELEI*(TRX*UVE-UVB)
     DUVE=DELEI*(TRX*DUVE-DUVB)
     IMAX=JRJ(JATOM)
     DO M=1,IMAX
        AE(M)=DELEI*(TRX*A(M)-AE(M))
        BE(M)=DELEI*(TRX*B(M)-BE(M))
     end DO
     !
     !     CALCULATE FUNCTION AT EI
     !
     CALL OUTWIN(REL,VR(1),R0(JATOM),DX(JATOM),JRJ(JATOM),EI,         &
          FL,UV,DUV,NODES,ZZ(jatom))
     CALL RINT13(A,B,A,B,OVLP,JATOM)
     TRX=1.0d0/SQRT(OVLP)
     P(l,1)=TRX*UV
     DP(l,1)=TRX*DUV
     IMAX=JRJ(JATOM)
     DO M=1,IMAX
        A(M)=TRX*A(M)
        B(M)=TRX*B(M)
     end DO
     !
     !     INSURE ORTHOGONALIZATION
     !
     CALL RINT13(A,B,AE,BE,CROSS,JATOM)
     TRY=-CROSS
     IMAX=JRJ(JATOM)
     DO M=1,IMAX
        AE(M)=(AE(M)+TRY*A(M))
        BE(M)=(BE(M)+TRY*B(M))
     end DO
     IMAX=JRJ(JATOM)
     DO I=1,IMAX
        RF1(I,l,1)=A(I)
        RF2(I,l,1)=B(I)
        RF1(I,l,2)=AE(I)
        RF2(I,l,2)=BE(I)
     end DO
     P(l,2)=UVE+TRY*P(l,1)
     DP(l,2)=DUVE+TRY*DP(l,1)
     CALL RINT13(AE,BE,AE,BE,PEI(l),JATOM)
     WRITE(unit_out,8) L,P(l,1),DP(l,1),P(l,2),DP(l,2)
  end DO lloop
!
! nun fur lo
!
  loloop: DO l=0,lomax
     irf=2
     iloloop: do jlo=1,ilo(l)
        if (lapw(l) .or. jlo/=1) then
           irf=irf+1
           DELE=2.0D-3
           DELEI=0.25D0/DELE
           FL=L
           EI=elo(l,jlo)/2.d0
           !
           !     CALCULATE FUNCTION AT EI
           IF(rlo(jlo,l)) THEN
              ei=elo(l,nloat)/2.d0
              kappa=l
              CALL diracout(rel,vr(1),r0(jatom),dx(jatom),jrj(jatom),    &
                   ei,kappa,uv,duv,nodes,zz(jatom))
              CALL dergl(a,b,r0(jatom),dx(jatom),jrj(jatom))
              DO m = 1, jrj(jatom)
                 r_m = r0(jatom)*exp(dx(jatom)*(m-1))
                 b(m) = b(m)*r_m/(2.d0*clight+(elo(l,jlo)- &
                      2.d0*vr(m)/r_m)/(2.d0*clight))
                 b(m)=b(m)*clight
              ENDDO
           ELSE
              CALL outwin(rel,vr(1),r0(jatom),dx(jatom),jrj(jatom),   &
                   ei,fl,uv,duv,nodes,zz(jatom))
           ENDIF

           CALL RINT13(A,B,A,B,OVLP,JATOM)
           TRX=1.0d0/SQRT(OVLP)
           P(l,irf)=TRX*UV
           DP(l,irf)=TRX*DUV
           IMAX=JRJ(JATOM)
           n_rad(l)=irf
           DO M=1,IMAX
              rf1(M,l,irf)=TRX*A(M)
              rf2(M,l,irf)=TRX*B(M)
           end DO

           CALL RINT13(rf1(1,l,1),rf2(1,l,1), &
                rf1(1,l,irf),rf2(1,l,irf),pi12lo,JATOM)
           CALL RINT13(rf1(1,l,2),rf2(1,l,2), &
                rf1(1,l,irf),rf2(1,l,irf),pe12lo,JATOM)
        end if
        call abc (l,jatom,pei(l),pi12lo,pe12lo,jlo,lapw(l))
     end do iloloop
  end DO loloop

  write(unit_out,651)(n_rad(l),l=0,lmax2)
651 format('number of rad. functions per L:',8I3)

  RETURN

5 FORMAT(10X,' ENERGY PARAMETERS ARE',7F7.2)
7 FORMAT(/10X,'ATOMIC PARAMETERS FOR ',A10/)
14 FORMAT(/11X,1HL,5X,4HU(R),10X, 5HU'(R),9X,5HDU/DE,8X,6HDU'/DE,6X,7HNORM-U')
8 FORMAT(10X,I2,5E14.6,5X,3I2)
1980 FORMAT(3X)
2000 FORMAT(15X,I3//)
2030 FORMAT(///)
2031 FORMAT(/)
END SUBROUTINE ATPAR


!!/---
!! Local Variables:
!! mode: f90
!! End:
!!\---
!!
!! Time-stamp: <2016-07-06 12:42:54 assman@faepop71.tu-graz.ac.at>
