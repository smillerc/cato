#!/bin/bash
set -e
set -u

source /software/intel/vtune_profiler/vtune-vars.sh intel64
export OMP_NUM_THREADS=6
cato_exe=/workspace/projects/cato/build/bin/cato.x

cp ${cato_exe} . &&
vtune --collect threading -- ./cato.x input.ini
