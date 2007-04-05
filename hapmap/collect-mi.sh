#!/bin/bash

# Collect calculated Mutual Information and print it as a table.

source config.sh

SUFFIX="$1"
SEP="||"
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
		MEAN_MI=$(cat ${POPULATION}${SUFFIX}/Chr22Results${STATES}States2/mean-coefficient-of-constraint.txt)
		echo -n $MEAN_MI | head -c 6
	done
	echo -n -e "$TRAIL_SEP"
	echo
done
