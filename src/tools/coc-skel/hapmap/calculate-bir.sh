#!/bin/bash

if [[ -r config.sh ]]
then
	source config.sh
else
	echo "Error: config.sh not found"
	exit 1
fi

function usage {
	echo "Usage: $0 [options]"
	echo "Options are:"
	echo "  --help       Print this message"
	echo "  --out <file> Output file name."
	echo "  --data <dir> Data directory, present directory by default."
	exit $1
}

# Default values
SOURCE_DATA_DIR="."
# Bir as in Bayesian Information Reward
OUT_FN="bir"

# Note that we use `"$@"' to let each command-line parameter expand to a 
# separate word. The quotes around `$@' are essential!
# We need TEMP as the `eval set --' would nuke the return value of getopt.
TEMP=$(getopt -o ha:sfc:o:m: --long help,data: \
     -n 'example.bash' -- "$@")

if [[ $? != 0 ]]; then echo "Terminating..." >&2 ; exit 1 ; fi

# Note the quotes around `$TEMP': they are essential!
eval set -- "$TEMP"

while true ; do
        case "$1" in
		-h|--help)
			usage 0 ;;

		-d|--data)
			SOURCE_DATA_DIR="$2"
			shift 2 ;;

		-o|--out)
			OUT_FN="$2"
			shift 2 ;;

		--) shift ; break ;;
                *) echo "Internal error! '$1'" ; exit 1 ;;
        esac
done
# echo "Remaining arguments:"
# for arg do echo '--> '"\`$arg'" ; done

if [[ -z "${OUT_FN}" ]]
then
	echo "No output file given."
	usage 1
fi


pushd "${SOURCE_DATA_DIR}"

for POPULATION in $POPULATIONS
do
	for STATES in $STATES_LIST
	do
		echo "${POPULATION}, ${STATES} states"
		R CMD BATCH --no-save --no-restore \
			--chromosome=Chr22         \
			--population=$POPULATION   \
			--states=$STATES           \
			bayesian-information-reward.R
		if [[ "$?" != "0" ]]
		then
			echo "R script has returned an error."
			tail "MutualInformation.Rout"
			exit 1
		fi
	done
done

popd
