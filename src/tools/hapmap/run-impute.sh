#!/bin/bash

POPULATIONS="Eur Afr Asian"
SRC_DATA_DIR="/home/maciej/chr22-data"
GENETIC_MAP="${SRC_DATA_DIR}/genetic_map_chr22.txt"

IMPUTE_DATA_BASE_DIR="$(pwd)"

function quiet {
"$@" > /dev/null
}

function run-impute {
for POPULATION in $POPULATIONS
do
	run_population $POPULATION
done
}

function display-and-run {
echo "$@"
"$@"
}

function run_population {
local POPULATION
POPULATION="$1"

IMPUTE_DATA_DIR="${IMPUTE_DATA_BASE_DIR}/${POPULATION}"
mkdir -p "${POPULATION}"

# Usage: impute-prepare.py [options]
# 
# Options:
#   -h, --help            show this help message and exit
#   -d DIR, --data-dir=DIR
#                         Data directory
#   -t FILE, --training=FILE
#                         Input training haplotypes file
#   -m FILE, --masked=FILE
#                         Input testing diplotypes file
#   -l FILE, --train-loci=FILE
#                         Input locus file
#   -a FILE, --haplo=FILE 
#                         Output haplo.txt file
#   -b FILE, --geno=FILE  Output geno.txt file
#   -c FILE, --legend=FILE
#                         Output legend.txt file

display-and-run ./impute-prepare.py --data-dir ${SRC_DATA_DIR}/${POPULATION}/chr22data \
	--haplo ${IMPUTE_DATA_DIR}/haplo.txt \
	--geno ${IMPUTE_DATA_DIR}/geno.txt \
	--legend ${IMPUTE_DATA_DIR}/legend.txt

# Run impute
display-and-run time impute \
	-h "${IMPUTE_DATA_DIR}/haplo.txt" \
	-l "${IMPUTE_DATA_DIR}/legend.txt" \
	-g "${IMPUTE_DATA_DIR}/geno.txt" \
	-m "${GENETIC_MAP}" \
	-Ne 11400
mv -v out "${IMPUTE_DATA_DIR}"/out.txt
mv -v info "${IMPUTE_DATA_DIR}"/info.txt
}

run-impute

