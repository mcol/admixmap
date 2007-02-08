#!/bin/bash

WORKING_DIR="legend"
OUT_FILE="rs_data.csv"

rm -f "$OUT_FILE"

for CHROMOSOME in $(seq 1 22)
do
	FILES="$WORKING_DIR/genotypes_chr${CHROMOSOME}_*.txt.gz"
	# echo $FILES
	for FILE in $FILES
	do
		echo -n "."
		CHROMOSOME=$(echo -n $FILE | cut -d_ -f2)
		zcat "$FILE" \
			| sed -e "s/.*/&\t$CHROMOSOME/" \
			>> "$OUT_FILE"
			# | awk '{ printf("%s\t'$CHROMOSOME'\n", $1);  }' \
	done
done
echo

echo "Sorting and removing duplicates it can take some time."
time \
cat "$OUT_FILE" | sort | uniq > tmp
mv tmp "$OUT_FILE"
