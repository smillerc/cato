#!/bin/bash

# Run this script inside a
export I_MPI_REMOVED_VAR_WARNING=0
export I_MPI_VAR_CHECK_SPELLING=0
export FOR_COARRAY_NUM_IMAGES=4

BUILD_DIR=../../../build

rm cato.x
rm -rf step*
cd ${BUILD_DIR} &&\
make -j &&\
cd - &&\
cp ${BUILD_DIR}/bin/cato.x . &&\
valgrind \
    --leak-check=full \
    --leak-resolution=high \
    --show-reachable=yes \
    --track-origins=yes \
    --log-file="log_mcheck.txt" ./cato.x input.ini
