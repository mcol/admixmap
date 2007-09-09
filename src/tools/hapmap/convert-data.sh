#!/bin/bash

# script to convert hapmixmap data used to assess prediction to IMPUTE format
# also writes a text file with commands for running impute, suitable for use with a taskfarm

WORKDIR=/ichec/work/ndlif006b
DATADIR=${WORKDIR}/data/chr22
HAPMIXMAP_DATADIR=${DATADIR}/hapmixmap
IMPUTE_DATADIR=${DATADIR}/impute

HAPMAP_MAPFILE=${DATADIR}/genetic_map_chr22.txt

PANELS="CEU YRI JPTCHB"
CONVERSION_TOOL=hapmix2impute

TASKFARM_FILE=impute-tasklist.txt
NE_VALUES="11418 17469 14269"

rm -f $TASKFARM_FILE

i=0
for PANEL in $PANELS
  do
  HAPMAP_LEGENDFILE=${WORKDIR}/genepi/hapmap/data/chr22/rawdata/${PANEL}/chr22_legend.txt
  HAPMIXMAP_TRAIN_GENOTYPES=${HAPMIXMAP_DATADIR}/train_genotypes.txt
  HAPMIXMAP_TEST_GENOTYPES=${HAPMIXMAP_DATADIR}/test_genotypes.txt
#  mkdir ${IMPUTE_DATADIR}/${PANEL}
#  ${CONVERSION_TOOL} -o=${IMPUTE_DATADIR}/${PANEL} -l=${HAPMAP_LEGENDFILE} \
#      -g=${HAPMIXMAP_TRAIN_GENOTYPES} -t=${HAPMIXMAP_TEST_GENOTYPES}\
#      -m=${HAPMAP_MAPFILE}
  echo "impute -h ${IMPUTE_DATADIR}/${PANEL}/haplo.txt -g ${IMPUTE_DATADIR}/${PANEL}/geno.txt -l ${IMPUTE_DATADIR}/${PANEL}/legend.txt ${IMPUTE_DATADIR}/${PANEL}/map.txt -o results/impute/${PANEL}_out.txt -i results/impute/${PANEL}/info.txt -Ne ${NE_VALUES}@$i " >>${TASKFARM_FILE}
i=$i+1
done