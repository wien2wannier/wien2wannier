          L_index=0
          index=0
          DO 119 L1=0,LMAX2
          DO 119 L2=0,LMAX2
          DO 119 LJ=0,LJMAX
           IF (MOD((L1+L2+LJ),2) .EQ. 1) GOTO 281
           IF ((L1+L2-LJ) .LT. 0) GOTO 281
           IF ((L1-L2+LJ) .LT. 0) GOTO 281
           IF ((-L1+L2+LJ) .LT. 0) GOTO 281
           L_index=L_index+1
           DO M1=-L1,L1
            DO MJ=-LJ,LJ                   ! compute overlap contribution of atom LATOM
             IF (ABS(M1+MJ).le.L2) THEN    ! overlap=conjg(alm(2))*alm(1)*gaunt(2,j,1)*rad_int(2,j,1)
              M2=M1+MJ
              index=index+1
              indexL(index)=L_index
              index1(index)=L1*(L1+1)+M1+1
              index2(index)=L2*(L2+1)+M2+1
              indexj(index)=LJ*(LJ+1)+MJ+1
              GNT(index)=GAUNT1(L2,LJ,L1,M2,MJ,M1)
              R_index=0
              DO irf1=1,nrf
              DO irf2=1,nrf
              R_index=R_index+1
              DO num1=1,NB
              DO num2=1,NB
              overlap(num2,num1,pair_index)=overlap(num2,num1,pair_index)+4*PI* &
              dconjg(ALM(INDEX2,num2,mu,irf2,2))*ALM(INDEX1,num1,mu,irf1,1)* &
              dconjg(YLB(indexj))*imag**LJ*PHSHEL* &
              ri_mat(R_index,L_index)*GAUNT1(L2,LJ,L1,M2,MJ,M1)
              ENDDO
              ENDDO
              ENDDO
              ENDDO
             ENDIF
            ENDDO
           ENDDO
 281       CONTINUE
 119      CONTINUE

