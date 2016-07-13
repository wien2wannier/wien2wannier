#!/bin/bash

### wien2wannier/test/PoC/SrVO3-PoC/SrVO3-PoC.sh
###
### Copyright 2016 Elias Assmann

set -e

w2wdir=$(realpath $(dirname $0)/../..)
x="$w2wdir/SRC/x_lapw -d"

me=SrVO3-PoC

exec >$(basename $PWD).log

init_lapw -b 2>&1

run_lapw 2>&1

$x lapw1 -band
lapw1 lapw1.def 2>&1

$w2wdir/SRC/write_insp_lapw -tmpl $WIENROOT/SRC_templates/case.insp

$x spaghetti
spaghetti spaghetti.def 2>&1

$w2wdir/SRC/prepare_w2wdir_lapw W

yes | clean_lapw

cd W 

$x kgen -fbz
echo -e '64\n0' | kgen kgen.def 2>&1

$w2wdir/SRC/write_inwf_lapw -bands 21 23 V:dt2g

$w2wdir/SRC_trig/write_win_backend <$w2wdir/SRC_templates/case.win \
    W.inwf W.struct W.klist W.klist_band > W.win

$w2wdir/SRC/wannier90_lapw -pp

$x lapw1
lapw1 lapw1.def 2>&1

$x w2w
$w2wdir/SRC_w2w/w2w w2w.def 2>&1

$w2wdir/SRC/wannier90_lapw

$w2wdir/SRC_trig/write_inwplot W <<EOF
0 0 0 1
1 0 0 1
0 1 0 1
0 0 1 1
10 10 10
EOF

$x wplot -wf 1
$w2wdir/SRC_wplot/wplot wplot.def

$w2wdir/SRC/wplot2xsf_lapw

## Time-stamp: <2016-07-13 16:51:27 assman@faepop71.tu-graz.ac.at>
