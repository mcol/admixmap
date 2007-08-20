#!/bin/bash
# 
# Generate shell scripts to run hapmixmap with Parkinson's data.

OUT_DIR="taskfarm-scripts"

function per_dir {
echo -n "${1}/${2}-chr${3}.sh"
}

function one_tf {
echo -n "${1}/taskfarm.sh"
}

# Strategy pattern (GoF), defining which function to call to get the
# file name.
#
# GET_FN="per_dir"
GET_FN="one_tf"

rm -f ${OUT_DIR}/*


for POP in ceu
do
	for CHR in $(seq 1 22)
	do
		# OUT_FILE="${OUT_DIR}/${POP}-chr${CHR}.sh"
		OUT_FILE="$($GET_FN ${OUT_DIR} ${POP} ${CHR})"
		( echo -n ". ${MODULESHOME}/init/sh ; module load pathscale ; "; \
		echo "./hmx-wrapper.pl --samples 1000 --burnin 80 --pop ${POP} --chromosome ${CHR}" ) >> "${OUT_FILE}"
		# echo "./hmx-wrapper.pl --samples 1000 --burnin 80 --pop ${POP} --chromosome ${CHR}" >> "${OUT_FILE}"
	done
done

