#!/bin/bash

# Error/unset flags
set -e
set -u

export I_MPI_REMOVED_VAR_WARNING=0
export I_MPI_VAR_CHECK_SPELLING=0
export FOR_COARRAY_NUM_IMAGES=1
export OMP_NUM_THREADS=6

cato_dir=../../../build
run_dir=`pwd`
rm -rf results_*

python generate_ic.py

if [ -f "cato.x" ]; then rm cato.x; fi
if [ -f "cato.error" ]; then rm cato.error; fi

# Get the executable
cd ${cato_dir} && make -j && cd ${run_dir} && cp ${cato_dir}/bin/cato.x .

echo "Creating both the symmetric and perturbed input files"
flux_solver='AUSM+-up_all_speed'
# flux_solver='FVLEG'

pert_ic='perturbed.h5'
pert_title='1D Perturbed'
pert_ini='perturbed.ini'

sym_ic='symmetric.h5'
sym_title='1D Symmetric'
sym_ini='symmetric.ini'

if [[ "$OSTYPE" == "linux-gnu" ]]; then
    sed_ver=sed
elif [[ "$OSTYPE" == "darwin"* ]]; then
    sed_ver=gsed
fi

echo "Setting the flux solver to: " $flux_solver
${sed_ver} -i '/^\[scheme\]$/,/^\[/ s/^.*flux_solver.*/flux_solver = "'"$flux_solver"'"/' _input.ini

# Remove existing
rm -f $pert_ini
rm -f $sym_ini

# Symmetric
echo "Making the symmetric input"
cp -v _input.ini ${sym_ini}
${sed_ver} -i '/^\[general\]$/,/^\[/ s/^.*title.*/title = "'"$sym_title"'"/' $sym_ini
${sed_ver} -i '/^\[initial_conditions\]$/,/^\[/ s/^.*initial_condition_file.*/initial_condition_file = "'"$sym_ic"'"/' $sym_ini

# Perturbed
echo "Making the perturbed input"
cp -v _input.ini ${pert_ini}
${sed_ver} -i '/^\[general\]$/,/^\[/ s/^.*title.*/title = "'"$pert_title"'"/' $pert_ini
${sed_ver} -i '/^\[initial_conditions\]$/,/^\[/ s/^.*initial_condition_file.*/initial_condition_file = "'"$pert_ic"'"/' $pert_ini

# Run
./cato.x $sym_ini 2>&1 | tee sym.out && mv results results_symmetric && \
./cato.x $pert_ini 2>&1 | tee pert.out && mv results results_perturbed
