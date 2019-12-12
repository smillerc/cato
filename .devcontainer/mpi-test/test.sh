#!/bin/sh
set -e

echo "--- Test MPI C installation ---"

printf "it should find mpicc... "
mpicc --version > /dev/null
echo ok

printf "it should find mpiexec... "
mpiexec --version > /dev/null
echo ok

printf "it should compile mpi_hello_world.c source... "
mpicc -o mpi_hello_world mpi_hello_world.c > /dev/null
echo ok

printf "it should run mpi_hello_world program successfully... "
mpirun -n `nproc` ./mpi_hello_world > /dev/null
echo ok

echo "--- Test MPI Fortran installation ---"

printf "it should find mpifort... "
mpifort --version > /dev/null
echo ok

printf "it should find mpiexec... "
mpiexec --version > /dev/null
echo ok

printf "it should compile mpi_hello_world.f90 source... "
mpifort -o mpi_hello_world_f08 mpi_hello_world.f90 > /dev/null
echo ok

printf "it should run mpi_hello_world program successfully... "
mpirun -n `nproc` ./mpi_hello_world_f08
echo ok
