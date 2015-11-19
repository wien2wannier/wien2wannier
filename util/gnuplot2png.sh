#!/bin/bash

if [ $# -ne 1 ]; then
	echo "Usage: $0 <wien2k seed>"
	echo "E.g.: $0 bands53-60"
	exit
fi 

SEED=$1

for FILE in *_x*.dat; do
BASENAME=`basename $FILE .dat`
gnuplot << EOT
set terminal png
set cbrange[0:2]
set xlabel "Z-Axis"
set ylabel "Y-Axis"
set output "${BASENAME}.png"
set pm3d map interpolate 10,10
splot "${FILE}
EOT
done
