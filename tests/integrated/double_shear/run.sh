#!/bin/bash
module purge
module load pfunit
module load cgns
module load intel/2019.3

export I_MPI_REMOVED_VAR_WARNING=0
export I_MPI_VAR_CHECK_SPELLING=0
export FOR_COARRAY_NUM_IMAGES=1

rm fvleg_2d.x
cd ../../../build &&\
make -j &&\
cd - &&\
cp ../../../build/bin/fvleg_2d.x . &&\
./fvleg_2d.x double_shear_input.ini
