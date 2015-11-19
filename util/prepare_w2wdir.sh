#!/bin/bash
#
#prepares a separate working directory for wien2wannier/wannier90
#N.Frohner, P.Wissgott 10/01/10
if [ "$1" == "-h" ]; then
  echo "prepares a separate working directory for wien2wannier/wannier90"
  echo "Usage: $0 [-c/-sp/-spc] case targetdir"
  echo "E.g.: $0 NiS2 k444"
  exit 0
fi

if [ $# -lt 2 ]; then
    echo "Usage: $0 [-c/-sp/-spc] case targetdir"
    echo "E.g.: $0 NiS2 k444"
    exit 0
fi 


if [ $# -ne 3 ]; then
   REQFILES="struct ksym scf2 output1 in1 vns vsp outputkgen dayfile klist_band spaghetti_ene"  
   SEED=$1
   TARGET=$2
else

  SEED=$2
  TARGET=$3
  if [ "$1" == "-c" ]; then
     REQFILES="struct ksym scf2 output1 in1c vns vsp outputkgen dayfile klist_band spaghetti_ene"
  elif [ "$1" == "-sp" ]; then
     REQFILES="struct ksym scf2up scf2dn  output1up output1dn in1 vnsup vnsdn vspup vspdn outputkgen dayfile klist_band spaghettiup_ene spaghettidn_ene"
  elif [ "$1" == "-spc" ]; then
     REQFILES="struct ksym scf2up scf2dn output1up output1dn in1c vnsup vnsdn vspup vspdn  outputkgen dayfile klist_band spaghettiup_ene spaghettidn_ene"
elif [ "$1" == "-spso" ]; then
     REQFILES="struct ksym scf2up scf2dn output1up output1dn in1c vnsup vnsdn vspup vspdn  outputkgen dayfile klist_band spaghettiup_ene spaghettidn_ene"
  else
    echo "Error: Unknown option"
    exit 0
  fi
fi




mkdir -p $TARGET
if [ ! -d $TARGET ]; then
		echo "Couldnt create target directory!"
		exit
fi

if [ ! -r $SEED.ksym ]; then
     cp $SEED.struct $SEED.ksym
     echo "Warning: copied $SEED.struct file to $SEED.ksym but did not change symmetry operations!"
      
fi

for FORMAT in $REQFILES; do
	REQFILE="${SEED}.${FORMAT}"
	if [ ! -r $REQFILE ]; then
		echo "Required File $REQFILE doesnt exist!"
		exit
	fi
	cp -i $REQFILE ${TARGET}/${TARGET}.${FORMAT}
        
done

if [ "$1" == "-c" ]; then
  write_w2wdef $TARGET
  mv w2w.def ${TARGET}/
elif [ "$1" == "-sp" ]; then
  write_w2wdef -up $TARGET
  mv w2w.def ${TARGET}/
elif [ "$1" == "-spc" ]; then
  write_w2wdef -up $TARGET
  mv w2w.def ${TARGET}/
else
  write_w2wdef $TARGET
  mv w2w.def ${TARGET}/
fi
