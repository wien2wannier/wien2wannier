#!/bin/bash
#testcase2 : SrVO3 spin-polarized
source get_seedname.sh

#Wien2k
cp $W2WROOT/templates/debug.struct $filename.struct
$WIENROOT/init -sp -b -numk 500
#cp ../templates/debug.klista debug.klist
runsp_lapw -cc 0.001

#bandstructure
cp $W2WROOT/templates/debug.klist_band $filename.klist_band
x lapw1 -up -band
x lapw1 -dn -band
cp $WIENROOT/SRC_templates/case.insp $filename.insp
update_FermiEnergy.sh -up
x spaghetti -up
update_FermiEnergy.sh -dn
x spaghetti -dn

#wien2wannier
cp $filename.struct $filename.ksym
cp $W2WROOT/templates/debug.klistb $filename.klist
cp $filename.klist $filename.klist_w90
x lapw1 -up
x lapw1 -dn
find_bands -up $filename -1 2
find_bands -dn $filename -1 2
cp $W2WROOT/templates/debug.w2win $filename.w2winup
cp $W2WROOT/templates/debug.w2win $filename.w2windn
cp $W2WROOT/templates/debug.outputkgen $filename.outputkgen
write_win -up $filename
write_win -dn $filename
wannier90.xsp -pp -up $filename
wannier90.xsp -pp -dn $filename
write_w2wdef -up $filename
w2wsp -up $filename
write_w2wdef -dn $filename
w2wsp -dn $filename

shift_energy -up $filename
shift_energy -dn $filename

#wannier90
wannier90.xsp -up $filename
wannier90.xsp -dn $filename

#wplot
write_wplotdef -up $filename
cp $W2WROOT/templates/debug.wplotin $filename.wplotinup
cp $W2WROOT/templates/debug.wplotin $filename.wplotindn
prepare_plotssp.sh -up $filename
prepare_plotssp.sh -dn $filename

#xcrysden
xsfAllsp.sh -up $filename
xsfAllsp.sh -dn $filename
xcrysden --xsf "$filename""_1.xsfup.gz" &
xcrysden --xsf "$filename""_1.xsfdn.gz" &

#rephase
rephase -up $filename
rephase -dn $filename



