#!/bin/bash
#cp debug.klist debug.klist_w90
cp $W2WROOT/templates/debug.klist_w2k debug.klist_w2k
cp debug.klist_w2k debug.klist
x lapw1
cp $W2WROOT/templates/debug.inop .
x optic
cp $W2WROOT/templates/debug.outputkgen_orig .
cp $W2WROOT/templates/debug.woptin debug.woptin
woptic -i 2 -inter debug >tmp.out
