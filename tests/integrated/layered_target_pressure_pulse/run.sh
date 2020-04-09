#!/bin/bash
#module purge
#module load pfunit
#module load cgns
#module load intel/2019.3

export I_MPI_REMOVED_VAR_WARNING=0
export I_MPI_VAR_CHECK_SPELLING=0
export FOR_COARRAY_NUM_IMAGES=1

python generate_ic.py
rm cato.x cato.error
#rm -rf step*
cd ../../../build &&\
make -j &&\
cd - &&\
cp ../../../build/bin/cato.x . &&\
./cato.x input.ini
echo
echo "Error(s) if any"
echo
cat cato.error
