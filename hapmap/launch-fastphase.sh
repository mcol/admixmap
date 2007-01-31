#!/bin/bash

FAST_PHASE="fastPHASE/fastphase.1.2.3.linux/fastPHASE"
FAST_PHASE="$(pwd)/$FAST_PHASE"
echo "fastPHASE binary: $FAST_PHASE"

source config.sh

STATES=8

# for POPULATION in $POPULATIONS
function populations {
for POPULATION in Eur
do
	# Convert the data.
	# Main data should be hiven as haploid with option -B, sourcing
	# from the mi_genotypes.txt and mi_loci.txt files, saved as
	# fastphase_haplotypes.inp. The main input file is
	# mi_cc_fastphase.inp
	SRC_DATA_DIR="$POPULATION/chr22data"
	WORKING_DATA_DIR="$POPULATION-fastPHASE/chr22data"
	WORKING_DIR="${POPULATION}-fastPHASE/Chr22Results${STATES}States2"
	mkdir -p "$WORKING_DIR" "$WORKING_DATA_DIR"
	# Copy the data
	for FILE in mi_genotypes.txt mi_loci.txt mi_cc_fastphase.inp mi_cc_observed_dput.txt mi_cc_index.txt
	do
		cp -v "$SRC_DATA_DIR/$FILE" "$WORKING_DATA_DIR"
	done
	# for FILE in mi_cc_observed_dput.txt
	# do
	# 	cp -v "$SRC_DATA_DIR" "$WORKING_DATA_DIR"
	# done
	pushd cut-loci
	perl convert-to-fastphase.pl \
		--haploid "../$WORKING_DATA_DIR/mi" \
		--fastphase "../$WORKING_DATA_DIR/fastphase_haplotypes.inp" \
		--no-loci-count
	popd
	pushd "$WORKING_DATA_DIR"
	$FAST_PHASE -T2 -C2 -K$STATES \
		-s50 \
		-bfastphase_haplotypes.inp \
		mi_cc_fastphase.inp
	popd
	# Extract the posterior distribution and save it as separate
	# file with values, in R `dget' format.
	python cut-loci/extract-posterior-probs.py \
		--sampled "$WORKING_DATA_DIR/fastphase_sampledHgivG.txt" \
		--output "$WORKING_DIR/PPGenotypeProbs.txt" \
		--args-file "$WORKING_DATA_DIR/mi_cc_index.txt" \
		-i 10
    echo "Running R for coefficient of constraint calculation."
	R CMD BATCH --no-save --no-restore \
		--chromosome=Chr22 \
		--population=${POPULATION}-fastPHASE \
		--states=$STATES \
		MutualInformation.R
done
}

populations
