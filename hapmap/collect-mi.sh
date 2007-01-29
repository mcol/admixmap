#!/bin/bash

# Collect calculated Mutual Information and print it as a table.

STATES_LIST=$(seq 2 2 18)
POPULATIONS="Eur Afr Asian"

echo -n "States"
for POPULATION in $POPULATIONS
do
	echo -n -e "\t"
	echo -n $POPULATION
done
echo

for STATES in $STATES_LIST
do
	echo -n "$STATES"
	for POPULATION in $POPULATIONS
	do
		echo -n -e "\t"
		MEAN_MI=$(cat $POPULATION/Chr22Results${STATES}States2/mean-coefficient-of-constraint.txt)
		echo -n $MEAN_MI
	done
	echo
done
