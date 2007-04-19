#!/bin/bash

. ${MODULESHOME}/init/sh
module load taskfarm

source config.sh

CHR22=chr22.pl

BASE_NAME="mutual-information"
TASK_FILE=${BASE_NAME}-tf.sh
PREPARATION_FILE=${BASE_NAME}-prep.sh
TIME_TOTAL="0"

rm -f "$TASK_FILE" "$PREPARATION_FILE"
touch "$TASK_FILE" "$PREPARATION_FILE"

# Predict execution time according to the regression model below
# Coefficients:
#                         Estimate Std. Error t value Pr(>|t|)    
# (Intercept)            2.960e+01  6.788e+00   4.361 8.83e-05 ***
# Samples:Haploid        2.473e-03  1.964e-04  12.596 1.67e-15 ***
# Samples:Haploid:States 5.111e-04  1.974e-05  25.891  < 2e-16 ***
# 
# New measurements (2007-04-10):
# (Intercept)           -6097.5046  1859.5035  -3.279 0.001253 ** 
# States                  595.8793   130.2158   4.576 8.86e-06 ***
# Samples                   0.6335     0.1850   3.425 0.000764 ***
# Haploid                  35.3335    12.0846   2.924 0.003907 ** 
predict_time() {
	P_SAM=$1
	P_HAP=$2
	P_STA=$3
	# old
	# EXPR="29.60 + $P_SAM * $P_HAP * 0.002473 + $P_SAM * $P_HAP * $P_STA * 0.0005111"
	# 2007-04-10
	EXPR="-6097.5046 + $P_STA * 595.8793 + $P_SAM * 0.6335 + $P_HAP + 35.3335"
	echo "$EXPR" | bc
}

for POPULATION in $POPULATIONS
do
	# rm -rfv "$POPULATION"
	DATADIR="$POPULATION/chr22data"
	echo >> "$PREPARATION_FILE" "mkdir -p -v \"$DATADIR\""
	echo >> "$PREPARATION_FILE" "cp -v ~/shared-genepi/hapmap/$DATADIR/phased_{genotypes,loci}.txt $DATADIR"
	echo -n >> "$PREPARATION_FILE" \
	perl $CHR22 \
		--mask-data \
		--pop=$POPULATION \
		--percent-indivs 17 \
		--percent-loci 30 \
		--limit-loci 5000 \
		--genotypes-file mi_genotypes.txt \
		--locus-file phased_loci.txt \
		--maskfile mi_cc_index.txt \
		--case-control-file mi_cc.txt
	echo >> "$PREPARATION_FILE"
	echo >> "$PREPARATION_FILE" 'if [ "$?" != "0" ]; then echo "Perl script returned an error."; exit 1; fi'
done

for POPULATION in $POPULATIONS
do
	# From the highest to the lowest number of states to minimize
	# cycle wastage when hitting the time limit.
	for STATES in $STATES_LIST
	do
		echo -n >> "$TASK_FILE" '. ${MODULESHOME}/init/sh'
		echo -n >> "$TASK_FILE" " ; "
		echo -n >> "$TASK_FILE" module load pathscale
		echo -n >> "$TASK_FILE" " ; "
		# Training, presumably a long run
		PREDICTED_TIME="$(predict_time 200 100 $STATES)"
		TIME_TOTAL="$TIME_TOTAL + $PREDICTED_TIME"
		# echo >> "$TASK_FILE" "# Predicted time: $PREDICTED_TIME"
		echo -n >> "$TASK_FILE" \
		perl $CHR22 \
			--pop $POPULATION \
			--states $STATES \
			--genotypes-file mi_genotypes.txt \
			--locus-file mi_loci.txt \
			--samples $SAMPLES \
			--burnin $(expr $SAMPLES / 5 ) \
			--re-run
		echo -n >> "$TASK_FILE" " ; "
		# FIXME: Number of individuals
		# predict_time $SAMPLES 100 $STATES
		PREDICTED_TIME=$(predict_time 200 100 $STATES)
		TIME_TOTAL="$TIME_TOTAL + $PREDICTED_TIME"
		# echo -n >> "$TASK_FILE" "# Predicted time: $PREDICTED_TIME"
		echo -n >> "$TASK_FILE" \
		perl $CHR22 \
			--pop $POPULATION \
			--states $STATES \
			--maskfile mi_cc_index.txt \
			--genotypes-file mi_merged_cc_train.txt \
			--locus-file mi_loci.txt \
			--samples ${TESTING_SAMPLES} \
			--burnin $(expr ${TESTING_SAMPLES} / 5 ) \
			--mutual-information \
			--re-run
		# new line
		echo >> "$TASK_FILE"
	done
done
echo -n "Total time: seconds: "
echo "($TIME_TOTAL)" | bc
echo -n "minutes: "
echo "($TIME_TOTAL) / 60.0" | bc
echo -n "hours: "
echo "($TIME_TOTAL) * 1.0 / 60.0 / 60.0" | bc
echo 
