#!/bin/bash

set -e
set -u

echo "Archiving the latest results into: $1"
mkdir $1
cp input.ini $1/input.ini
cp generate_ic.py $1/generate_ic.py
cp energy_input.dat $1/energy_input.dat
cp step* $1
