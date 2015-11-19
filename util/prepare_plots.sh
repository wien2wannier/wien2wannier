#!/bin/bash
#plots all Wannier orbitals of a series 
#P.Wissgott 10/01/10
if [ $# -ne 1 ]; then
    echo "Usage: prepare_plots.sh case"
    echo "e.g.: prepare_plots.sh SrVO3"
    exit 0
elif [ "$1" == "-h" ]; then
    echo "plots all Wannier orbitals of a series"
    echo "the number of Wannier orbitals is read from wout file"
    echo "the squared values |w_m(r)|^2 are stored to case_m.psink"
    echo "the phases phi(w_m(r)) are stored to case_m.psiarg"
    echo "Usage: prepare_plots.sh case"
    echo "e.g.: prepare_plots.sh SrVO3"
    exit 0
fi

filename=$1
filename_wout="${filename}"".wout"
echo "read number of wannier functions from $filename_wout"
tmp=$(grep "Number of Wannier Functions" $filename_wout)
echo $tmp
num_wann=$(echo ${tmp:47:30})
echo "number of Wannier functions is $num_wann"
for (( idx = 1;  idx <= num_wann;idx++    )) do
   
  echo "compute Wannier function $idx ......."
  wplot ${filename} $idx
  cp "${filename}.psiarg" "${filename}""_""$idx"".psiarg"
  cp "${filename}"".psink" "${filename}""_""$idx"".psink"
done
