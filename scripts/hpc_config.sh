#!/bin/bash

module purge
module load intel/2019.5
module load hdf5/1.8.17/b2
module rm intel/17U4 # part of the hdf5 mod
module load gcc/9.1.0

which ifort
which mpiifort
rm -rf CMake* CTest* Testing bin include lib src tests cmake* Make*
find ../src  -name "*.mod" -delete
CC=gcc FC=mpiifort cmake .. -DCMAKE_BUILD_TYPE="Release"

make -j
