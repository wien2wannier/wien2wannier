#!/bin/bash

### wien2wannier/test/PoC/PoC.sh
###
### Copyright 2016 Elias Assmann

set -e

cd $1

exec >$1.log

init_lapw -b $(cat $1.opt-init)    2>&1
run_lapw -p -ec 0.1                2>&1

x lapw1 -band -p                   2>&1
write_insp_lapw
x spaghetti -p                     2>&1

yes | clean_lapw

prepare_w2wdir_lapw W
cd W

init_w2w -b $(cat ../$1.opt-w2w)

x lapw1 -p                         2>&1
x w2w -p                           2>&1
x wannier90                        2>&1

if [[ -s ../$1.inwplot ]]; then
    cp ../$1.inwplot W.inwplot
else
    write_inwplot W <<EOF
0 0 0 1
1 0 0 1
0 1 0 1
0 0 1 1
10 10 10
EOF
fi

x wplot -wf 1 -p                   2>&1
wplot2xsf_lapw                     2>&1

