#!/bin/bash

export CUDA_VISIBLE_DEVICES=0

MODEL_NAME="1a99"
#min
pmemd -O -i min-imp.in -o "$MODEL_NAME"_min.out -c "$MODEL_NAME".crd -p "$MODEL_NAME".top -r "$MODEL_NAME"_min.rst -ref "$MODEL_NAME".crd
#heat1
pmemd.cuda -O -i heat-imp.in -o "$MODEL_NAME"_heat1.out -c "$MODEL_NAME"_min.rst -p "$MODEL_NAME".top -r "$MODEL_NAME"_heat1.rst -x "$MODEL_NAME"_heat1.mdcrd -ref "$MODEL_NAME"_min.rst 
#equilibrate 1
pmemd.cuda -O -i density-imp.in -o "$MODEL_NAME"_equil1.out -c "$MODEL_NAME"_heat1.rst -p "$MODEL_NAME".top -r "$MODEL_NAME"_equil1.rst -x "$MODEL_NAME"_equil1.mdcrd -ref "$MODEL_NAME"_heat1.rst
#equilibrate 2
pmemd.cuda -O -i equil-imp.in -o "$MODEL_NAME"_equil2.out -c "$MODEL_NAME"_equil1.rst -p "$MODEL_NAME".top -r "$MODEL_NAME"_equil2.rst -x "$MODEL_NAME"_equil2.mdcrd -ref "$MODEL_NAME"_equil1.rst

#md
pmemd.cuda -O -i md-imp.in -o "$MODEL_NAME"_md1.out -c "$MODEL_NAME"_equil2.rst -p "$MODEL_NAME".top -r "$MODEL_NAME"_md1.rst -x "$MODEL_NAME"_md1.mdcrd
