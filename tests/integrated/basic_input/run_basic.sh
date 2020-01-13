#!/bin/bash

rm cato.x
cd ../../../build &&\
make -j &&\
cd - &&\
cp ../../../build/bin/cato.x . &&\
./cato.x basic.ini
