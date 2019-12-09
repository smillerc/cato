#!/bin/bash

rm fvleg_2d.x
cd ../../../build &&\
make -j &&\
cd - &&\
cp ../../../build/bin/fvleg_2d.x . &&\
./fvleg_2d.x double_shear_input.ini
