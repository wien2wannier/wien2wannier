#!/bin/bash
seed=$1
savename=$2
cp $seed.w2win $savename.w2win
cp $seed.win $savename.win
cp $seed.mmn $savename.mmn
cp $seed.amn $savename.amn
cp $seed.nnkp $savename.nnkp
cp $seed.klist_w90 $savename.klist_w90
cp $seed.ham $savename.ham
cp $seed"_hr.dat" $savename"_hr.dat"
cp $seed"_band.dat" $savename"_band.dat"
cp $seed.chk $savename.chk
cp $seed.eig $savename.eig
cp $seed.wout $savename.wout