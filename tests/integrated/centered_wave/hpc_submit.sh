#!/bin/bash

# Set up variables
JOBNAME=${PWD##*/}
JOBDIR=`pwd`
EXEC="/b1/smil/cato/build/bin/cato.x"
QUEUE=s1
NODES=1
CPUS=12
NP=$(($NODES*$CPUS))

# The lines below shouldn't need to be modified to often

### Make the script
cat > pbs_submit.sh << EOF
#!/bin/bash
#PBS -l select=${NODES}:ncpus=${CPUS}:mpiprocs=${CPUS}
#PBS -q ${QUEUE}
#PBS -N ${JOBNAME:0:8}
#PBS -j oe
#PBS -m a
#PBS -m e

cd ${JOBDIR}
cp ${EXEC} .

module purge
module load intel/2019.5
module load hdf5/1.8.17/b2
module rm intel/17U4 # part of the hdf5 mod
module load gcc/9.1.0

./cato.x input.ini
exit
EOF

### Submit the job
qsub pbs_submit.sh

