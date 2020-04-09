#!/bin/bash

export I_MPI_REMOVED_VAR_WARNING=0
export I_MPI_VAR_CHECK_SPELLING=0
export FOR_COARRAY_NUM_IMAGES=1

if [[ "$OSTYPE" == "linux-gnu" ]]; then
    cato_dir=$(readlink -f ../../../build)
elif [[ "$OSTYPE" == "darwin"* ]]; then
    cato_dir=$(greadlink -f ../../../build)
fi

run_dir=`pwd`

python generate_ic.py

if [ -f "cato.x" ]; then rm cato.x; fi
if [ -f "cato.error" ]; then rm cato.error; fi
rm -rf step*

cd ${cato_dir} && make -j && \
    cd ${run_dir} && \
    cp ${cato_dir}/bin/cato.x . &&\
    ./cato.x input.ini

if [ -f "cato.error" ]; then
    echo
    echo "Error(s) if any"
    echo
    cat cato.error
fi
