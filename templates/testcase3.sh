#!/bin/bash
#testcase2 : SrVO3 spin-polarized incl SO
source get_seedname.sh

#Wien2k
cp $W2WROOT/templates/debug2.struct $filename.struct
cp $W2WROOT/templates/debug2.inso $filename.inso
cp $W2WROOT/templates/debug2.in1 $filename.in1
cp $W2WROOT/templates/debug2.in2c $filename.in2c
$WIENROOT/init -sp -b -numk 500
runsp_lapw -cc 0.001 -so

#bandstructure
cp $W2WROOT/templates/debug.klist_band $filename.klist_band
x lapw1 -up -band
x lapw1 -dn -band
x lapwso -up -band
cp $WIENROOT/SRC_templates/case.insp $filename.insp
update_FermiEnergy.sh -up
x spaghetti -up -so
update_FermiEnergy.sh -dn
x spaghetti -dn -so

#wien2wannier
cp $filename.struct $filename.ksym
cp $W2WROOT/templates/debug.klistb $filename.klist
cp $filename.klist $filename.klist_w90
x lapw1 -up
x lapw1 -dn
x lapwso -up
find_bands -soup $filename -1 2
find_bands -sodn $filename -1 2
cp $W2WROOT/templates/debug2.w2win $filename.w2winup
cp $W2WROOT/templates/debug2.w2win $filename.w2windn
cp $W2WROOT/templates/debug.outputkgen $filename.outputkgen
write_win -up $filename
wannier90.xsp -pp -up $filename
w2wsp -so $filename

shift_energy -up $filename
cp $filename.eigup $filename.eigso

#wannier90
wannier90.xso $filename

#wplot
write_wplotdef -soup $filename
cp $W2WROOT/templates/debug.wplotin $filename.wplotinup
prepare_plotsso.sh -up $filename
write_wplotdef -sodn $filename
cp $W2WROOT/templates/debug.wplotin $filename.wplotindn
prepare_plotsso.sh -dn $filename

#xcrysden
xsfAllsp.sh -soup $filename
xsfAllsp.sh -sodn $filename
xcrysden --xsf "$filename""_1.xsfup.gz" &
xcrysden --xsf "$filename""_1.xsfdn.gz" &




