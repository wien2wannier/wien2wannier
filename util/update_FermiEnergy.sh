#!/bin/bash

### wien2wannier/util/update_FermiEnergy.sh
###
###    Writes Fermi energy from scf2 file to insp file
###
### Copyright 2009-2012 Philipp Wissgott
### Copyright      2013 Elias Assmann

if [ $# -eq 1 ] && [ "$1" == "-h" ]; then
   echo "writes Fermi energy from scf2 file to insp file"
   echo "Usage: update_FermiEnergy.sh [-up/dn]"
   exit 0
fi

case=$(basename $PWD)
insp=$case.insp

if [ $# -eq 0 ]; then
    scf=$case.scf2
elif [ $# -eq 1 ] && [ "$1" == "-up" ]; then
    scf=$case.scf2up
elif [ $# -eq 1 ] && [ "$1" == "-dn" ]; then
    scf=$case.scf2dn
else
   echo "Unknown option: $1"
   exit 0
fi

echo "looking in $scf for Fermi energy"
EF=$(grep ":FER" $scf | tail -n1 |
     perl -ne 'print m/.*=\s*([0-9.]+)/')

echo "Fermi energy is $EF Ry"

if [[ ! -s $insp ]]
then
    cp $WIENROOT/SRC_templates/case.insp $insp
fi

perl -i~ -pe \
'$.==9 or next;
m/^(\s*\S*\s+)(\S+)(\s*)(.*)$/;
$_ = $1 . '$EF' .
     q( ) x (length($2) + length($3) - '${#EF}') .
     $4 . qq(\n)' $insp
