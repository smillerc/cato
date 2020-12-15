#!/bin/bash

# Error checking
set -u

export FOR_COARRAY_NUM_IMAGES=1

cato_dir=../../../build
run_dir=`pwd`
rm -rf results
python generate_ic.py

if [ -f "cato.x" ]; then rm cato.x; fi
if [ -f "std.err" ]; then rm std.err; fi

cd ${cato_dir} && make -j && \
    cd ${run_dir} && \
    cp ${cato_dir}/bin/cato.x . &&\
    cafrun -np ${FOR_COARRAY_NUM_IMAGES} ./cato.x input.ini

if [ -f "std.err" ]; then
    echo
    echo "Error(s) if any"
    echo
    cat std.err
fi

#python view_results.py
