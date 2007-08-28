#!/bin/bash

PRIORS="
irw-2-arp-1.2-0.05-1-0.5-afpp-0.25-1 \
irw-2-arp-1.2-0.05-1-0.5-afpp-2-8 \
irw-2-arp-1.2-0.05-2-afpp-0.25-1 \
irw-2-arp-1.2-0.05-2-afpp-2-8 \
irw-2-arp-12-0.5-1-0.5-afpp-0.25-1 \
irw-2-arp-12-0.5-1-0.5-afpp-2-8 \
irw-2-arp-12-0.5-2-afpp-0.25-1 \
irw-2-arp-12-0.5-2-afpp-2-8 \
"


WORKDIR=/ichec/work/ndlif006b/maciej
PREFIX=irw-2
RESULTSDIR=Chr22Results6States2
POPULATIONS="Eur Afr Asian"
SEEDS=1

for POPULATION in $POPULATIONS
do
  for PRIOR in $PRIORS
    do
    for SEED in $SEEDS
      do
#      echo "${POPULATION}, seed ${SEED}"
      "R --vanilla --args  \
datadir=${WORKDIR}/chr22-data/${POPULATION}/chr22data \
resultsdir=${WORKDIR}/${PREFIX}-${PRIOR}-${SEED}/hapmap/${POPULATION}/${RESULTSDIR} \
<../hapmap/bayesian-information-reward.R >BIR.Rout" 

      if [[ "$?" != "0" ]]
	  then
	  echo "R script has returned an error."
	  tail "BIR.Rout"
	  exit 1
      fi
    done
  done
done

