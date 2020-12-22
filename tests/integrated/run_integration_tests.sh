#!/bin/bash

# Run integration tests
export I_MPI_REMOVED_VAR_WARNING=0
export I_MPI_VAR_CHECK_SPELLING=0
export FOR_COARRAY_NUM_IMAGES=6
export OMP_NUM_THREADS=1
export UCX_LOG_LEVEL=error

set -u
set -e

# Set the vendor to gnu or intel
vendor=$1
cato_dir=../../../build
cato_exe=${cato_dir}/bin/cato.x

if [ $vendor == "intel" ]; then
    run=''
else
    run="cafrun -np ${FOR_COARRAY_NUM_IMAGES}"
fi

# 1D Tests (serial mode)
for test in sod_1d shu_osher_shocktube
do
    echo "Test: " ${test}
    cd ${test} && \
    ${cato_exe} input.ini && \
    python view_results.py && \
    cd ..
done

# 2D Tests (use coarrays)
for test in shock_noise sedov implosion
do
    echo "Test: " ${test}
    cd ${test} && \
    ${run} ${cato_exe} input.ini && \
    python view_results.py && \
    cd ..
done