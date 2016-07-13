#!/bin/bash

### wien2wannier/test/PoC/SrVO3-PoC/SrVO3-PoC.sh
###
### Copyright 2016 Elias Assmann

set -e

exec >$(basename $PWD).log

init_lapw -b   2>&1
run_lapw       2>&1

x lapw1 -band  2>&1
write_insp_lapw
x spaghetti    2>&1

yes | clean_lapw

prepare_w2wdir_lapw W
cd W 

init_w2w -bands 21 23 -proj V:dt2g

x lapw1        2>&1
x w2w          2>&1
x wannier90    2>&1

write_inwplot W <<EOF
0 0 0 1
1 0 0 1
0 1 0 1
0 0 1 1
10 10 10
EOF

x wplot -wf 1  2>&1
wplot2xsf_lapw 2>&1

## Time-stamp: <2016-07-13 18:31:14 assman@faepop71.tu-graz.ac.at>
