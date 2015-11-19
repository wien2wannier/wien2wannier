#!/bin/bash
#testcase1 : SrVO3 simple calculation
source get_seedname.sh

#Wien2k
cp $W2WROOT/templates/debug.struct $filename.struct
$WIENROOT/init -b -numk 500
#cp ../templates/debug.klista debug.klist
run_lapw -cc 0.001

#bandstructure
cp $W2WROOT/templates/debug.klist_band $filename.klist_band
x lapw1 -band
cp $WIENROOT/SRC_templates/case.insp $filename.insp
update_FermiEnergy.sh
x spaghetti

#wien2wannier
cp $filename.struct $filename.ksym
cp $W2WROOT/templates/debug.klistb $filename.klist
cp $filename.klist $filename.klist_w90
x lapw1
find_bands $filename -1 1
cp $W2WROOT/templates/debug.w2win $filename.w2win
cp $W2WROOT/templates/debug.outputkgen $filename.outputkgen
write_win $filename
wannier90.x -pp $filename
write_w2wdef $filename
w2w $filename

shift_energy $filename

#wannier90
wannier90.x $filename

#wplot
write_wplotdef $filename
cp $W2WROOT/templates/debug.wplotin $filename.wplotin
prepare_plots.sh $filename

#xcrysden
xsfAll.sh $filename
xcrysden --xsf "$filename""_1.xsf.gz" &

#rephase
rephase -w $filename





