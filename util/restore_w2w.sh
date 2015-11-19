#!/bin/bash
seed=$1
savename=$2
cp $savename.w2win $seed.w2win
cp $savename.win $seed.win
cp $savename.mmn $seed.mmn
cp $savename.amn $seed.amn
cp $savename.nnkp $seed.nnkp
cp $savename.klist_w90 $seed.klist_w90
cp $savename.ham $seed.ham
cp $savename"_hr.dat" $seed"_hr.dat"
cp $savename"_band.dat" $seed"_band.dat"
cp $savename.chk $seed.chk
cp $savename.eig $seed.eig
cp $savename.wout $seed.wout