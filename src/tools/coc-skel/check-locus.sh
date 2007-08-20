#!/bin/bash
# 
# Print out summary for a column of hapmixmap genotypes data file.
# Usage:
# bash check-locus.sh <file> <snp-id>
#
# The script prints all the values found in a certain column in
# hapmixmap data file, along with the snp identifier.

GT_FILE="$1"
SNP_ID="$2"

function summarize_column {
	GT_FILE1="$1"
	COL_NO1="$2"
	echo "COL_NO: ${COL_NO1}"
	awk "{ print \$${COL_NO1};}" < "${GT_FILE1}" | sort | uniq -c
}

COL_NO="$(head -n 1 ${GT_FILE} | head -n 1 | tr '\t' '\n' | grep -n ${SNP_ID} |  cut -d: -f1)"
summarize_column "${GT_FILE}" "$(( ${COL_NO} - 1 ))" "${SNP_ID}"
summarize_column "${GT_FILE}" "$(( ${COL_NO}     ))" "${SNP_ID}"
summarize_column "${GT_FILE}" "$(( ${COL_NO} + 1 ))" "${SNP_ID}"
