!!! wien2wannier/SRC_w2w/l2mmn.f

SUBROUTINE l2MMN(NB,num_kpts,NNTOT,LJMAX)
  USE param
  USE struct 
  USE xa
  USE xa3
  USE bessel
  USE amn_mmn
  USE pairs
  use lolog, only: nlo,nlov,nlon,loor,ilo,lapw,n_rad
  IMPLICIT REAL(R8) (A-H,O-Z)
  INTEGER          pair_index, r_index
  REAL(R8)           KX1,KY1,KZ1                       ! k-points in u.c. coordinates k and k+b
  COMPLEX(C16)       YLB((LMAX2+1)*(LMAX2+1))       ! spherical harmonics expansion of b
  COMPLEX(C16)       PHSHEL,tmp,tmp1
  COMMON /GENER/   BR1(3,3),BR2(3,3)              ! transformation between u.c. and cartesian coordinates
  COMMON /ATSPDT/  P(0:LMAX2,nrf),DP(0:LMAX2,nrf) ! radial function and its slope at RMT
  COMMON /RADFU/   RF1(NRAD,0:LMAX2,nrf),RF2(NRAD,0:LMAX2,nrf)  ! radial functions large and small component
  common /loabc/   alo(0:lomax,nloat,nrf)         ! use for local orbitals

  !...............................
  complex(C16), allocatable :: alm(:,:,:,:),blm(:,:,:,:)
  !...............................

  integer :: k1_prog_itvl
  k1_prog_itvl = min(max(num_kpts/10, 1), 100)

  !------------------------------------------------------------------     

  TWOPI=2.D0*PI  
  FOURPI=4*PI

  OVERLAP = 0

  !     ------------------
  !     LOOP FOR ALL ATOMS
  !     ------------------

  READ(unit_vsp,2032) ISCF    
  LFIRST=1
  atoms: DO JATOM=1,NAT 
     write(unit_out, "(/, '===== atom', I5, ' /', I5, ' =====' /)") jatom, nat

     ALLOCATE(ALM(NB,NRF,(LMAX2+1)*(LMAX2+1),MULT(JATOM)),  &
          BLM(NB,NRF,(LMAX2+1)*(LMAX2+1),MULT(JATOM)))
     talm=0.
     tmeas1=0.
     tmeas2=0.
     call cputim(t1)
     IF(JATOM.GT.1) LFIRST=LFIRST+MULT(JATOM-1)                      
     ITAP=30+JATOM                                                     
     itape=unit_vector
     jtape=unit_vsp
     rewind(itape)
     CALL ATPAR(JATOM,LFIRST,itape,jtape)      ! calculate radial functions for atoms JATOM
     FAC=4.0D0*PI*RMT(JATOM)**2/SQRT(VOL)
     rewind(itape)

     !................................
     pair_index=0
     k1loop: DO k1=1,num_kpts      !   loop over k1
        pair_index=pair_index+1
        kkk=KP(pair_index)
        KX1=XK(kkk)
        KY1=YK(kkk)
        KZ1=ZK(kkk)
        !	write(92,*)'k-point:',kkk,KX1,KY1,KZ1
        CALL almgen(ALM,JATOM,LFIRST,NB,KKK)
        k2loop: DO k2=1,NNTOT
           if (k2.gt.1) pair_index=pair_index+1
           kkk=KPB(pair_index)
           call cputim(tt0)
           CALL almgen(BLM,JATOM,LFIRST,NB,KKK)
           BLM=conjg(BLM)
           call cputim(tt1)
           BX=XK(kkk)-KX1+BQX(pair_index)                 ! calculate b=k2-k1 add BQ if going around BZ
           BY=YK(kkk)-KY1+BQY(pair_index)
           BZ=ZK(kkk)-KZ1+BQZ(pair_index)
           !	 write(92,*)BX,BY,BZ
           BK(1)=BX*BR1(1,1)+BY*BR1(1,2)+BZ*BR1(1,3)      ! transform to cartesian coordinates
           BK(2)=BX*BR1(2,1)+BY*BR1(2,2)+BZ*BR1(2,3)
           BK(3)=BX*BR1(3,1)+BY*BR1(3,2)+BZ*BR1(3,3)
           BM=SQRT(BK(1)*BK(1)+BK(2)*BK(2)+BK(3)*BK(3))
           CALL RADINT(JATOM,LJMAX,BM)                  ! compute radial intergrals <R(r)|j_(|b|*r)|R'(r)>
           CALL YLM (BK,LJMAX,YLB)                      ! computer Y_lm(b)
           indexj=0
           DO LJ=0,LJMAX
              DO MJ=-LJ,LJ
                 indexj=indexj+1
                 YLB(indexj)=conjg(YLB(indexj))*(0,1)**LJ
              ENDDO
           ENDDO
           call cputim(tt2)
           muloop: DO mu=1,MULT(JATOM)
              latom=lfirst-1+mu
              ARG1=BX*POS(1,LATOM)*TWOPI
              ARG2=BY*POS(2,LATOM)*TWOPI
              ARG3=BZ*POS(3,LATOM)*TWOPI
              PHSHEL=EXP((0,1)*(ARG1+ARG2+ARG3))*FOURPI 

              ! overlap=conjg(alm(2))*alm(1)*gaunt(2,j,1)*rad_int(2,j,1)
              L_index=0
              l1loop: DO L1=0,LMAX2
                 l2loop: do L2=0,LMAX2
                    ljloop: DO LJ = abs(L1-L2), min(L1+L2, LJMAX), 2
           IF (MOD((L1+L2+LJ),2) .EQ. 1) cycle
           IF ((L1+L2-LJ).LT.0.OR.(L1-L2+LJ).LT.0.OR.(-L1+L2+LJ).LT.0) cycle
                       L_index=L_index+1
                       m1loop: DO M1=-L1,L1
                          mjloop: DO MJ = max(-LJ, -L2-M1), min(LJ, L2-M1)
                             M2=M1+MJ                   ! abs(m2) <= l2 !
                             index1=L1*(L1+1)+M1+1
                             index2=L2*(L2+1)+M2+1
                             indexj=LJ*(LJ+1)+MJ+1
                             tmp=YLB(indexj)*PHSHEL*GAUNT1(L2,LJ,L1,M2,MJ,M1)
                             R_index=0
                             DO irf1=1,n_rad(L1)
                                do irf2=1,n_rad(L2)
                                   R_index=R_index+1
                                   tmp1=ri_mat(R_index,L_index)*tmp
                                   do num1=1,NB
                                      do num2=1,NB
                                         overlap(num2,num1,pair_index) = &
                                              overlap(num2,num1,pair_index) + &
                                              BLM(num2,irf2,INDEX2,mu) * &
                                              ALM(num1,irf1,INDEX1,mu) * &
                                              tmp1
                                      end do
                                   end do
                                ENDDO
                             ENDDO
                          ENDDO mjloop
                       ENDDO m1loop
                    ENDDO ljloop
                 ENDDO l2loop
              ENDDO l1loop
           ENDDO muloop
           call cputim(tt3)
           talm=talm+tt1-tt0
           tmeas1=tmeas1+tt2-tt1
           tmeas2=tmeas2+tt3-tt2
        END DO k2loop

        if (mod(k1, k1_prog_itvl) == 0) &
             write(unit_out, "('k1=', I5, ' /', I5, ' (', I3, '%)')") &
             &    k1, num_kpts, (100*k1)/num_kpts
     END DO k1loop
     DEALLOCATE(ALM,BLM)
     call cputim(t2)
     talm=talm*(NNTOT+1)/NNTOT
     write(unit_out,*) 'CPU time used for atom ',JATOM,' ='&
          ,t2-t1,talm,tmeas1,tmeas2

        if (mod(k1, k1_prog_itvl) == 0) &
             write(unit_out, "('k1=', I5, ' /', I5, ' (', I3, '%)')") &
             &    k1, num_kpts, (100*k1)/num_kpts
  END DO atoms

  ! ....END LOOP OVER ALL ATOMS     

  RETURN

2032 FORMAT(50X,I2,//)                                                        
   contains
     subroutine ptime(descr, unit)
       character(len=*), intent(in), optional :: descr
       integer,          intent(in), optional :: unit

       character(len=*), parameter :: fmt = "('Times for ', A, T33, '(sec):', &
            & F8.3, ' wall;', F9.3)"!, ' cpu =', F8.3, ' user +', F8.3, ' sys')"

       real(r8),     save :: cputime1, cputime2
       integer,    save :: walltime1, walltime2, count_rate
       integer,    save :: default_lun
       integer          :: lun

       if (.not. present(descr)) then
          call cpu_time(cputime1)
          call system_clock(walltime1, count_rate)

          if (present(unit)) default_lun=unit

          return
       end if

       if (present(unit)) then
          lun=unit
       else
          lun=default_lun
       end if

       call cpu_time(cputime2)
       call system_clock(walltime2)

       write(lun, fmt) descr, real(walltime2-walltime1)/count_rate, &
            & (cputime2-cputime1)

       walltime1 = walltime2
       cputime1  = cputime2
     end subroutine ptime
END SUBROUTINE l2MMN


!!/---
!! Local Variables:
!! mode: f90
!! End:
!!\---
!!
!! Time-stamp: <2015-12-28 15:13:21 assman@faepop36.tu-graz.ac.at>
