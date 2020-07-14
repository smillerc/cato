#!/bin/bash

# Error/unset flags
set -e
set -u

export I_MPI_REMOVED_VAR_WARNING=0
export I_MPI_VAR_CHECK_SPELLING=0
export FOR_COARRAY_NUM_IMAGES=1
export OMP_NUM_THREADS=6

cato_dir=../../../build
base_dir=`pwd`
rm -rf results*

if [ -f "cato.x" ]; then rm cato.x; fi
# Get the executable
cd ${cato_dir} && make -j && cd ${base_dir} && cp ${cato_dir}/bin/cato.x .

# Select the version of sed to use
if [[ "$OSTYPE" == "linux-gnu" ]]; then
    sed_ver=sed
elif [[ "$OSTYPE" == "darwin"* ]]; then
    sed_ver=gsed
fi

for solver in {'AUSM+-up_all_speed','FVLEG'}; do
    for pert in {'symmetric','perturbed'}; do
        folder=${solver}_${pert}
        rm -rf $folder
        mkdir -p $folder

        python generate_ic.py $pert
        mv initial_conditions.h5 $folder/initial_conditions.h5
        cp -v pressure_input.dat $folder/pressure_input.dat

        # Symmetric
        echo "Running the ${folder} case"
        cd ${folder}
        cp -v ../_input.ini input.ini

        ${sed_ver} -i '/^\[scheme\]$/,/^\[/ s/^.*flux_solver.*/flux_solver = "'"$solver"'"/' input.ini
        ${sed_ver} -i '/^\[general\]$/,/^\[/ s/^.*title.*/title = "'"$folder"'"/' input.ini
        ../cato.x input.ini 2>&1 | tee log.out
        cd ..
    done
done

echo "Done!"
