#!/bin/bash

export GMON_OUT_PREFIX 'gmon.out'

BUILD_DIR=../../../build

rm cato.x
rm -rf step*
cd ${BUILD_DIR} &&\
make -j &&\
cd - &&\
cp ${BUILD_DIR}/bin/cato.x . && ./cato.x input.ini

# Profile
gprof -A cato.x gmon.out > profile_results.txt
