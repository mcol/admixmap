changequote([,])dnl
define([PBS_SCRIPT], [#!/bin/bash
#PBS -N $2
#PBS -l nodes=4:ppn=2
#PBS -l walltime=00:05:00
#PBS -A ndlif006b

PBS_O_WORKDIR=$1/hapmap

. ${MODULESHOME}/init/sh
module load taskfarm

echo Working directory is ${PBS_O_WORKDIR}
cd "${PBS_O_WORKDIR}"

# Run a taskfarm
mpiexec taskfarm ./mutual-information-tf.sh])dnl
