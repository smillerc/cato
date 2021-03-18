#!/bin/bash

# Error/unset flags
set -e
set -u

export I_MPI_REMOVED_VAR_WARNING=0
export I_MPI_VAR_CHECK_SPELLING=0
export FOR_COARRAY_NUM_IMAGES=3
export OMP_NUM_THREADS=1

cato_dir=${HOME}/cato/build
PY=python

# Select the version of sed to use
if [[ "$OSTYPE" == "linux-gnu" ]]; then
    sed_ver=sed
    elif [[ "$OSTYPE" == "darwin"* ]]; then
    sed_ver=gsed
fi

for pert in {'symmetric','perturbed'}; do
    folder=${pert}
    rm -rf ${folder}
    mkdir -p ${folder}
    
    echo "Running the ${folder} case"
    cd ${folder}
    run_dir=`pwd`
    cp -v ../*.csv .
    cp -v ../input.ini input.ini
    cp -v ../generate_ic.py .
    ${sed_ver} -i '/^\[general\]$/,/^\[/ s/^.*title.*/title = "'"$folder"'"/' input.ini
    
    ${PY} generate_ic.py ${pert}
    
    cd ${cato_dir} && make -j && \
    cd ${run_dir} && \
    cp ${cato_dir}/bin/cato.x . &&\
    cafrun -np ${FOR_COARRAY_NUM_IMAGES} ./cato.x input.ini
    cd ..
done

echo "Done!"
