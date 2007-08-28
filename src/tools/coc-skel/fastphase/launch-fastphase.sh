#!/bin/bash

FAST_PHASE="fastPHASE/fastphase.1.2.3.linux/fastPHASE"
FAST_PHASE="$(pwd)/${FAST_PHASE}"
echo "fastPHASE binary: ${FAST_PHASE}"

source config.sh

STATES=8

function populations {
for POPULATION in $POPULATIONS
# for POPULATION in Asian
do
	# Convert the data.
	# Main data should be hiven as haploid with option -B, sourcing
	# from the mi_genotypes.txt and mi_loci.txt files, saved as
	# fastphase_haplotypes.inp. The main input file is
	# mi_cc_fastphase.inp
	SRC_DATA_DIR="${POPULATION}/chr22data"
	WORKING_DATA_DIR="${POPULATION}-fastPHASE/chr22data"
	WORKING_DIR="${POPULATION}-fastPHASE/Chr22Results${STATES}States2"
	mkdir -p "${WORKING_DIR}" "${WORKING_DATA_DIR}"
	# Copy the data
	for FILE in mi_genotypes.txt mi_loci.txt mi_cc_fastphase.inp mi_cc_observed_dput.txt mi_cc_index.txt
	do
		cp -v "$SRC_DATA_DIR/$FILE" "$WORKING_DATA_DIR"
	done
	# for FILE in mi_cc_observed_dput.txt
	# do
	# 	cp -v "$SRC_DATA_DIR" "$WORKING_DATA_DIR"
	# done
	perl convert-to-fastphase.pl \
		--haploid "../${WORKING_DATA_DIR}/mi" \
		--fastphase "../${WORKING_DATA_DIR}/fastphase_haplotypes.inp" \
		--no-loci-count
	pushd "$WORKING_DATA_DIR"
	# -M2 options has the same effect as `fixed mixture proportions'
	# in hapmixmap.
	# 
	# If this option is used, data file has to contain the P-line,
	# containing the relative loci positions.
	#
	# CAUTION: if the P-line is not present in the data, fastPHASE
	# will silently ignore the -M2 option.
	FP_OPTIONS="-T20 -C50 -K${STATES} \
		-s1000 \
		-p \
		-M2 \
		-bfastphase_haplotypes.inp \
		mi_cc_fastphase.inp"
	echo > fastphase_cli_options.txt ${FP_OPTIONS}
	${FAST_PHASE} ${FP_OPTIONS}
	if [ ! "$?" = "0" ]
	then
		echo "fastPHASE returned an error (dir $(pwd))."
		exit 1
	fi
	popd

	# Number of individuals needs to be hardcoded, as fastPHASE
	# doesn't output it in the file with sampled genotypes.
	if [ "${POPULATION}" = "Asian" ]
	then
		NO_INDIVS=15
	else
		NO_INDIVS=10
	fi

	# Extract the posterior distribution and save it as separate
	# file with values, in R `dget' format.
	python extract-posterior-probs.py \
		--sampled "${WORKING_DATA_DIR}/fastphase_sampledHgivG.txt" \
		--output "${WORKING_DIR}/PPGenotypeProbs.txt" \
		--args-file "${WORKING_DATA_DIR}/mi_cc_index.txt" \
		-i ${NO_INDIVS}
	# This file occupies a LOT of space. We better compress it
	# nicely, in the background.
	rm -v "${WORKING_DATA_DIR}/fastphase_sampledHgivG.txt.bz2"
	nice bzip2 "${WORKING_DATA_DIR}/fastphase_sampledHgivG.txt" &
	if [ "$?" != "0" ]
	then
		echo "Python script returned an error."
		exit 1
	fi
	echo "Running R for coefficient of constraint calculation."
	R CMD BATCH --no-save --no-restore \
		--chromosome=Chr22 \
		--population=${POPULATION}-fastPHASE \
		--states=${STATES} \
		../hapmap/MutualInformation.R
done
}

populations
