#!/bin/bash

module purge
#module load opencoarrays-2.7.1-gcc-9.2.0-usk2akd
module load hdf5-1.10.5-gcc-9.2.0-rvcszud
export HDF5_ROOT='/Users/smil/opt/spack/opt/spack/darwin-mojave-broadwell/gcc-9.2.0/hdf5-1.10.5-ipjrokunwzzib3fjbedkpgtjcsudujn6'
# export PFUNIT_DIR="/Users/smil/opt/PFUNIT-4.0"
export PFUNIT_DIR="/Users/smil/opt/pfunit_serial/PFUNIT-4.0"


rm -rf CMake* CTest* Testing bin include lib src tests cmake* Make*
find ../src  -name "*.mod" -delete
CC=gcc-9 FC=gfortran-9 cmake .. -DCMAKE_BUILD_TYPE="Release"
