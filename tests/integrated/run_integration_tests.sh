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

if [[ "$OSTYPE" == "linux-gnu" ]]; then
    results_dir=`readlink -f ./test_results`
elif [[ "$OSTYPE" == "darwin"* ]]; then
    results_dir=`greadlink -f ./test_results`
fi

rm -rf ${results_dir} && mkdir -p ${results_dir}

if [ $vendor == "intel" ]; then
    run=''
else
    run="cafrun -np ${FOR_COARRAY_NUM_IMAGES}"
fi

# 1D Tests (serial mode)
for test in sod_1d shu_osher_shocktube
do
    echo "Test: " ${test}
    cd ${test} && rm -rf results \
    ${cato_exe} input.ini && \
    python view_results.py && cp -v *.png ${results_dir} && cd ..
done

# 2D Tests (use coarrays)
for test in shock_noise sedov kelvin_helmholtz implosion
do
    echo "Test: " ${test}
    cd ${test} && rm -rf results \
    ${run} ${cato_exe} input.ini && \
    python view_results.py && cp -v *.png ${results_dir} && cd ..
done