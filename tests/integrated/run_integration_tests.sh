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
ls ../..

if [[ "$OSTYPE" == "linux-gnu" ]]; then
    results_dir=`readlink -f ./test_results`
    cato_dir=`readlink -f ../../build`
elif [[ "$OSTYPE" == "darwin"* ]]; then
    results_dir=`greadlink -f ./test_results`
    cato_dir=`greadlink -f ../../build`
fi

cato_exe=${cato_dir}/bin/cato.x
echo "CATO exe: ${cato_exe}"
rm -rf ${results_dir} && mkdir -p ${results_dir}

if [ $vendor == "intel" ]; then
    run=''
else
    run="cafrun -np ${FOR_COARRAY_NUM_IMAGES}"
fi

# 1D Tests (serial mode)
for test in sod_1d shu_osher_shocktube
do
    cd ${test} && rm -rf results && \
    echo "Test: " ${test} `pwd` && ls && \
    python generate_ic.py && \
    ${cato_exe} input.ini && \
    echo "Post-processing" && \
    python view_results.py && cp -v *.png ${results_dir} && cd ..
done

# 2D Tests (use coarrays)
for test in shock_noise sedov kelvin_helmholtz implosion
do
    cd ${test} && rm -rf results && \
    echo "Test: " ${test} `pwd` && ls && \
    python generate_ic.py && \
    ${run} ${cato_exe} input.ini && \
    echo "Post-processing" && \
    python view_results.py && cp -v *.png ${results_dir} && cd ..
done