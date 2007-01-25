#!/bin/bash

. ${MODULESHOME}/init/sh
module load taskfarm

CHR22=chr22.pl

TASKFILE=mutual-information-tf.sh
POPULATIONS="Eur Afr Asian"


rm -f "$TASKFILE"
touch "$TASKFILE"
SAMPLES=500

for POPULATION in $POPULATIONS
do
	# rm -rfv "$POPULATION"
	DATADIR="$POPULATION/chr22data"
	mkdir -p -v "$DATADIR"
	cp -v ~/shared-genepi/hapmap/$DATADIR/phased_{genotypes,loci}.txt $DATADIR
	echo -n >> "$TASKFILE" \
	perl $CHR22 \
		--mask-data \
		--percent-indivs 17 \
		--percent-loci 30 \
		--limit-loci 5000 \
		--genotypes-file mi_genotypes.txt \
		--locus-file phased_loci.txt \
		--maskfile mi_cc_index.txt \
		--case-control-file mi_cc.txt
	echo >> "$TASKFILE"
done

for POPULATION in $POPULATIONS
do
	for STATES in 2
	do
		# Training, presumably long run
		echo -n >> "$TASKFILE" \
		perl $CHR22 \
			--pop $POPULATION \
			--states $STATES \
			--genotypes-file mi_genotypes.txt \
			--locus-file mi_loci.txt \
			--samples $SAMPLES \
			--burnin $(expr $SAMPLES / 5 ) \
			--re-run
		echo -n >> "$TASKFILE" " ; "
		echo -n >> "$TASKFILE" \
		perl $CHR22 \
			--pop $POPULATION \
			--states $STATES \
			--maskfile mi_cc_index.txt \
			--genotypes-file mi_merged_cc_train.txt \
			--locus-file mi_loci.txt \
			--samples 50 \
			--burnin 10 \
			--mutual-information \
			--re-run
		echo >> "$TASKFILE"
	done
done
