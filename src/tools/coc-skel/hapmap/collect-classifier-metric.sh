#!/bin/bash

# Collect calculated metric and print it as a table.
#

# METRIC="coefficient-of-constraint"
METRIC="bayesian-information-reward"

source config.sh

SUFFIX="$1"
SEP="$2"
INT_SEP=" $SEP "
HEAD_SEP="$SEP "
TRAIL_SEP=" $SEP"

echo -n -e "$HEAD_SEP"
echo -n "States"
for POPULATION in $POPULATIONS
do
	echo -n -e "$INT_SEP"
	echo -n $POPULATION
done
echo -n -e "$TRAIL_SEP"
echo

for STATES in $STATES_LIST
do
	echo -n -e "$HEAD_SEP"
	echo -n "$STATES"
	for POPULATION in $POPULATIONS
	do
		echo -n -e "$INT_SEP"
		MEAN_MI_FILE_NAME="${POPULATION}${SUFFIX}/Chr22Results${STATES}States2/mean-${METRIC}.txt"
		if [ -r "${MEAN_MI_FILE_NAME}" ]
		then
			MEAN_MI="$(cat ${MEAN_MI_FILE_NAME})"
		else
			MEAN_MI="NA"
		fi
		echo -n ${MEAN_MI} | head -c 6
	done
	echo -n -e "$TRAIL_SEP"
	echo
done
