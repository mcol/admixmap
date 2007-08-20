#!/bin/bash

HAPMAP_DIR="/home/maciej/ucd/genedc-bzr/src/hapmap/download"
HAPMIXMAP_DIR="hapmap-hapmixmap"

function launch_python {
	POP=$1
	CHR=$2
	PADDED_CHR="$(printf "%02d" ${CHR})"
	python hapmap2hapmixmap.py \
		-a "$HAPMAP_DIR/genotypes-chr${CHR}-${POP}-r21-nr-fwd-phased" \
		-b "$HAPMAP_DIR/genotypes-chr${CHR}-${POP}-r21-nr-fwd-legend.txt" \
		-s "$HAPMAP_DIR/genotypes-chr${CHR}-${POP}-r21-nr-fwd-sample.txt" \
		-c "${HAPMIXMAP_DIR}/${POP}-chr${PADDED_CHR}-gt.txt" \
		-d "${HAPMIXMAP_DIR}/${POP}-chr${PADDED_CHR}-loci.txt"
}

# Executing two conversions in parallel. Processing chromosome 1 and
# 2 at the same time isn't a good idea because of the memory
# requirements.  Here, chromosome 1 is processed along with chromosome
# 22, 2 with 21 and so on. This leads to some in-efficiency (only one
# processor working), but at least doesn't lead to memory exhaustion
# (which happens with a 2GB machine).
# for POP in ceu jpt-chb yri
for POP in yri
do
	for CHR in $(seq 1 11)
	do
		echo "hapmap2hapmixmap.py: ${POP}, chr${CHR}"
		launch_python ${POP} ${CHR} &
		PID1="$1"
		launch_python ${POP} $(( 22 - ${CHR} + 1 )) &
		PID2="$1"
		wait $PID1
		wait $PID2
	done
done
