#!/bin/bash

PRIORS="
arp-1.2-0.05-1-0.5-afpp-0.25-1 \
arp-1.2-0.05-1-0.5-afpp-2-8 \
arp-1.2-0.05-2-afpp-0.25-1 \
arp-1.2-0.05-2-afpp-2-8 \
arp-12-0.5-1-0.5-afpp-0.25-1 \
arp-12-0.5-1-0.5-afpp-2-8 \
arp-12-0.5-2-afpp-0.25-1 \
arp-12-0.5-2-afpp-2-8 \
"

METRIC="bayesian-information-reward"
WORKDIR=/ichec/work/ndlif006b/maciej
PREFIX=irw-2
RESULTSDIR=Chr22Results6States2
POPULATIONS="Eur Afr Asian"
SEEDS=1
SEP="|"

echo "Population${SEP}Prior${SEP}Seed${SEP}${BIR}"
for POPULATION in $POPULATIONS
do
  for PRIOR in $PRIORS
    do
    for SEED in $SEEDS
      do
#      echo "${POPULATION}, seed ${SEED}"
      "R --vanilla --args  \
datadir=${WORKDIR}/chr22-data/${POPULATION}/chr22data \
resultsdir=${WORKDIR}/${PREFIX}-${PRIOR}-s${SEED}/hapmap/${POPULATION}/${RESULTSDIR} \
<../hapmap/bayesian-information-reward.R >BIR.Rout" 

      if [[ "$?" != "0" ]]
	  then
	  echo "R script has returned an error."
	  tail "BIR.Rout"
      fi

      echo -n -e "$POPULATION$SEP$PRIOR$SEP$SEED$SEP"
      MEAN_MI_FILE_NAME="${WORKDIR}/${PREFIX}-${PRIOR}-s${SEED}/hapmap/${POPULATION}/${RESULTSDIR}/mean-${METRIC}.txt"
      if [ -r "${MEAN_MI_FILE_NAME}" ]
	  then
	  MEAN_MI="$(cat ${MEAN_MI_FILE_NAME})"
      else
	  MEAN_MI="NA"
      fi
      echo -n ${MEAN_MI} | head -c 6
    done
  done
done

