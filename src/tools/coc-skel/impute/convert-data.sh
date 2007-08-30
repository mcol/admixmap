#!/bin/bash

#script to convert hapmixmap data used to assess prediction to IMPUTE format
# on Walton

WORKDIR=/ichec/work/ndlif006b
DATADIR=$(WORKDIR)/maciej/chr22-data
PANELS="Eur Afr Asian"
CONVERSION_TOOL=hapmix2impute

for PANEL in PANELS
do
HAPMIXMAP_DATADIR=${DATADIR}/${PANEL}/chr22data
IMPUTE_DATADIR=${DATADIR}/${PANEL}/impute-chr22data
HAPMAP_LEGENDFILE=${WORKDIR}/genepi/hapmap/${PANEL}/chr22_legend.txt
HAPMAP_MAPFILE=${WORKDIR}/genepi/hapmap/${PANEL}/genetic_map_chr22.txt
HAPMIXMAP_TRAIN_GENOTYPES=${HAPMIXMAP_DATADIR}/
HAPMIXMAP_TEST_GENOTYPES=${HAPMIXMAP_DATADIR}/
mkdir ${IMPUTE_DATADIR}
${CONVERSION_TOOL} -o=${IMPUTE_DATADIR} -l=${HAPMAP_LEGENDFILE} \
-g=${HAPMIXMAP_TRAIN_GENOTYPES} -t=${HAPMIXMAP_TEST_GENOTYPES}\
-m=${HAPMAP_MAPFILE}
done