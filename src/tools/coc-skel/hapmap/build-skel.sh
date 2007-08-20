#!/bin/bash

SRC_DIR=../genepi-old

for FILE in calculate-mi.sh chr22.pl collect-mi.sh config.sh count-fields.pl \
	cut-loci FastPhaseConverter.R generate-taskfarm-file.sh launch-fastphase.sh \
	launch-taskfarm.pbs log-time.sh maskGenotypesFunctions.R \
	mutual-information-prep.sh MutualInformation.R mutual-information-tf.sh
do
	cp -av ${SRC_DIR}/hapmap/${FILE} hapmap
done
