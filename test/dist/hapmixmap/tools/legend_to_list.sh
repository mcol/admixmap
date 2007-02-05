#!/bin/bash

WORKING_DIR="legend"
OUT_FILE="rs_data.csv"

rm -f "$OUT_FILE"

for CHROMOSOME in $(seq 1 22)
do
	FILES="$WORKING_DIR/genotypes_chr${CHROMOSOME}_*.txt.gz"
	echo $FILES
	for FILE in $FILES
	do
		CHROMOSOME=$(echo -n $FILE | cut -d_ -f2)
		zcat "$FILE" \
			| awk '{ printf("%s\t'$CHROMOSOME'\n", $1);  }' \
			>> "$OUT_FILE"
	done
done
