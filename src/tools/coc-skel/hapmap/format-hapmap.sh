#!/bin/bash

# Convert hapmap source data into hapmixmap format, retaining the split
# by one chromosome

POP=ceu
CHR=15

HAPMAP_SRC_DIR="hapmap-source"
DST_DIR="hapmap-hapmixmap"

for CHR in $(seq 22)
do
	for POP in ceu yri jpt-chb
	do
		HAPMIXMAP_G="${DST_DIR}/chr${CHR}-${POP}-genotypes.txt"
		echo -n "${HAPMIXMAP_G}..."

		echo -n -e "ID\t" > "${HAPMIXMAP_G}"

		cat "${HAPMAP_SRC_DIR}/genotypes-chr${CHR}-${POP}-r21-nr-fwd-legend.txt" \
			| grep -v position \
			| awk '{ print $1; }' \
			| tr '\n' '\t' | sed -e 's/\t$//' >> "${HAPMIXMAP_G}"
		echo >> "${HAPMIXMAP_G}"

		cat "${HAPMAP_SRC_DIR}/genotypes-chr${CHR}-${POP}-r21-nr-fwd-sample.txt" \
			| cut -d" " -f1 | nl > tmp1.txt

		cat "${HAPMAP_SRC_DIR}/genotypes-chr${CHR}-${POP}-r21-nr-fwd-phased" \
			| nl | tr 01 12 > tmp2.txt

		join tmp1.txt tmp2.txt | sed -e 's/^[[:digit:]]\+\s\+//g' >> "${HAPMIXMAP_G}"

		rm tmp1.txt tmp2.txt
		echo "done."
	done
done
