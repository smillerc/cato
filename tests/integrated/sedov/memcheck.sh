#!/bin/bash
#module purge
#module load pfunit
#module load cgns
#module load intel/2019.3

export I_MPI_REMOVED_VAR_WARNING=0
export I_MPI_VAR_CHECK_SPELLING=0
export FOR_COARRAY_NUM_IMAGES=1

rm cato.x
rm -rf step*
cd ../../../build &&\
make -j &&\
cd - &&\
cp ../../../build/bin/cato.x . &&\
valgrind --leak-check=full --leak-resolution=high --show-reachable=yes --track-origins=yes --log-file="log_mcheck.txt" ./cato.x input.ini
