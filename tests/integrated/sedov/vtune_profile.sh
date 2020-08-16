#!/bin/bash
set -e
set -u

source /software/intel/vtune_profiler/vtune-vars.sh intel64

THREADS=6
export OMP_DISPLAY_ENV=TRUE
export OMP_NUM_THREADS=${THREADS}
export KMP_AFFINITY=verbose,granularity=fine,compact,0,0

cato_dir=../../../build
cato_exe=${cato_dir}/bin/cato.x

cp ${cato_exe} . &&
vtune --collect threading -- ./cato.x input.ini
