#!/bin/bash

set -u
set -e

for test in sod_1d shock_noise sedov double_shear riemann_2d_prob_1
do
    echo "Test: " ${test}
    cd ${test} && ls && cd ..

done