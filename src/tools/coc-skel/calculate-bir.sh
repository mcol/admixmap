#!/bin/bash

source=$1
#arg should be a file containing a definition of DIRS where you want to calculate BIR
#and optionally PREFIX and SUFFIX

WORKDIR=/ichec/work/ndlif006b
PREFIX=irw-2
RESULTSDIR=Chr22Results6States2
POPULATIONS="Eur Afr Asian"
SEEDS=1

for POPULATION in $POPULATIONS
  do
  for DIR in $DIRS
    do
    echo "${DIR}"
    "R --vanilla --args  \
datadir=${WORKDIR}/maciej/chr22-data/${POPULATION}/chr22data \
resultsdir=${WORKDIR}/${PREFIX} \
<../hapmap/bayesian-information-reward.R >BIR.Rout" 

if [[ "$?" != "0" ]]
    then
    echo "R script has returned an error."
    tail "BIR.Rout"
    exit 1
fi

  done
done

