#!/bin/bash

#script to convert hapmixmap data used to assess prediction to IMPUTE format

WORKDIR=/ichec/work/ndlif006b
DATADIR=$(WORKDIR)/data/chr22
HAPMIXMAP_DATADIR=${DATADIR}/hapmixmap
IMPUTE_DATADIR=${DATADIR}/impute

HAPMAP_MAPFILE=${DATADIR}/genetic_map_chr22.txt

PANELS="Eur Afr Asian"
CONVERSION_TOOL=hapmix2impute

for PANEL in PANELS
  do
  HAPMAP_LEGENDFILE=${WORKDIR}/genepi/hapmap/data/chr22/rawdata/${PANEL}/chr22_legend.txt
  HAPMIXMAP_TRAIN_GENOTYPES=${HAPMIXMAP_DATADIR}/train_genotypes.txt
  HAPMIXMAP_TEST_GENOTYPES=${HAPMIXMAP_DATADIR}/test_genotypes.txt
  mkdir ${IMPUTE_DATADIR}/${PANEL}
  ${CONVERSION_TOOL} -o=${IMPUTE_DATADIR}/${PANEL} -l=${HAPMAP_LEGENDFILE} \
      -g=${HAPMIXMAP_TRAIN_GENOTYPES} -t=${HAPMIXMAP_TEST_GENOTYPES}\
      -m=${HAPMAP_MAPFILE}
done