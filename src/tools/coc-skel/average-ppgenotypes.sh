#!/bin/bash -x
# -*- coding: UTF-8 -*-

source hapmap/config.sh

function usage() {
echo
echo "Usage: bash average-ppgenotypes.sh --output <dir> <src1> [ <src2> ... ]"
echo "Options:"
echo "    --help             This message         "
echo "    --output <dir>     Destination directory"
echo "    --purge            Purge the destination directory"
echo
exit "$1"
}

TEMP=$(getopt -o ho:p --long help,output:,purge \
     -n 'example.bash' -- "$@")
if [[ $? != 0 ]]; then echo "Terminating..." >&2 ; exit 1 ; fi
eval set -- "$TEMP"

while true ; do
        case "$1" in
		-h|--help)
			usage 0
			shift
			;;
		-p|--purge)
			PURGE="1"
			shift
			;;
		-o|--output)
			OUT_DIR="$2"
			shift 2 ;;
		--) shift ; break ;;
                *) echo "Internal error! '$1'" ; exit 1 ;;
        esac
done

if [[ -z "${OUT_DIR}" ]]
then
	usage 1
fi

if [[ "${PURGE}" = "1" && -d "${OUT_DIR}" ]]
then
	echo "Removing existing \`${OUT_DIR}' directory"
	rm -rf "${OUT_DIR}"
fi

if [[ -d "${OUT_DIR}" ]]
then
	echo "Directory '${OUT_DIR}' exists already. Delete it first by hand."
	exit 1
fi

echo "Creating the output directory: \`${OUT_DIR}'"
mkdir -p "${OUT_DIR}"/hapmap/{Eur,Asian,Afr}/Chr22Results8States2

tar cf - hapmap/config.sh hapmap/calculate-mi.sh hapmap/collect-mi.sh \
	hapmap/*.R \
	| tar xvfp - -C "${OUT_DIR}"

tar cf - -C "$1" hapmap/{Asian,Eur,Afr}/chr22data/mi_cc_observed_dput.txt \
	| tar xvfp - -C "${OUT_DIR}"

# Assemble an R script to average over all GP's found

for POPULATION in $POPULATIONS
do
	echo "Population ${POPULATION}"
	echo > tmp.R "# vim:set ft=r:"

	N=0
	GP_LIST=""
	for arg
	do
		echo '--> '"\`$arg'"
		N=$(( $N + 1 ))
		echo >> tmp.R "GP${N} = dget(\"${arg}/hapmap/${POPULATION}/Chr22Results8States2/PPGenotypeProbs.txt\")"
		GP_LIST="$(echo $GP_LIST GP${N})"
	done

	GP_LIST="$(echo ${GP_LIST} | sed -e 's/ / + /g')"

	echo "Number of dirs to average: ${N}"
	echo >> tmp.R "GP = (" ${GP_LIST} ") / ${N}"
	echo >> tmp.R "dput(GP, \"${OUT_DIR}/hapmap/${POPULATION}/Chr22Results8States2/PPGenotypeProbs.txt\")"
	R CMD BATCH --no-save --no-restore tmp.R
	if [[ "$?" != "0" ]]
	then
		echo "R script failed."
		cat tmp.Rout
		exit 1
	fi
done

pushd "${OUT_DIR}/hapmap"
bash calculate-mi.sh
bash collect-mi.sh "" '||'
popd
# tree "${OUT_DIR}"
