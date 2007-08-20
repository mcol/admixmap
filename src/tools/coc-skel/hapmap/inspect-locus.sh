#!/bin/bash

function usage {
echo "Usage: $0 <file> <snp_id>"
echo
exit $1
}

function col_no {
FILE_NAME="$1"
SNP_ID="$2"
cat "${FILE_NAME}" | head -n 1 | tr '\t' '\n' | grep -n ^${SNP_ID}$ | cut -d: -f1
}

function inspect_by_col_no {
FILE_NAME="$1"
COL_NO1="$2"
echo "COL_NO1: '${COL_NO1}'"
cat  "${FILE_NAME}" | tr '\t' ' ' | cut -d " " -f ${COL_NO1} | sort | uniq -c
}

SNP_ID="$2"
FILE_NAME="$1"

[[ ! -z "${FILE_NAME}" ]] || usage 1
[[ ! -z "${SNP_ID}" ]] || usage 1

COL_NO="$(col_no ${FILE_NAME} ${SNP_ID})"

inspect_by_col_no "${FILE_NAME}" $(( ${COL_NO} - 1))
inspect_by_col_no "${FILE_NAME}" $(( ${COL_NO} ))
inspect_by_col_no "${FILE_NAME}" $(( ${COL_NO} + 1))
