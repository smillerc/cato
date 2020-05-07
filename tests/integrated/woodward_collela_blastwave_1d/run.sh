#!/bin/bash
#module purge
#module load pfunit
#module load cgns
#module load intel/2019.3

export I_MPI_REMOVED_VAR_WARNING=0
export I_MPI_VAR_CHECK_SPELLING=0
export FOR_COARRAY_NUM_IMAGES=1
export OMP_NUM_THREADS=2

cato_dir=../../../build
run_dir=`pwd`
rm -rf results
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
