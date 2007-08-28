#!/bin/bash
# -*- coding: UTF-8 -*-
#
# Set up a set of analyses
#
# This file is part of Genepi, genetic data analysis software.
# Copyright (C) 2007 Maciej Bliziński
# 
# Genepi is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# 
# Genepi is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with Genepi; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

########################################################################
# Script configuration
########################################################################

ISO_DATE="$(date +%Y-%m-%d)"
TIMESTAMP="$(date +%Y%m%d-%H%M)"

MAIN_DIRECTORY="shared-genepi/maciej"
RESULTS_DIR="html"
HAPMIXMAP_RES_SUBDIR="Chr22Results6States2"

PURGE_DIR=yes

PBS_FILE="launch-taskfarm.pbs"

RESDIRS="Eur Afr Asian"

# SOURCE_DATA_DIR="${HOME}/genepi-old/hapmap"
# SOURCE_DATA_DIR="$HOME/shared-genepi/maciej/chr22-tuned-2/hapmap"
SOURCE_DATA_DIR="" # no default

# Additional configuration:
# pbs.m4
# Please note that this script submits only 3 tasks per job. It submits
# the minimal number of nodes: 4, which is twice too much, but that's
# the minimum.

########################################################################
# Script body
########################################################################

function usage() {
echo "Usage: bash manage-analyses.sh [ options ]"
echo "Options:"
echo "    --action <action>         Action, see below"
echo "    --config <config>         Configuration file"
echo "    --data <dir>              Source data directory (no default)"
echo "  [ --force                ]  Delete directories when they exist."
echo "                              This can destroy existing results!"
echo "  [ --submit               ]  Submit tasks to cluster"
echo "  [ --output <dir>         ]  HTML output directory (default: \`html')"
echo "  [ --main-directory <dir> ]  Main working directory (default: \`shared-genepi/maciej')"
echo ""
echo "Action are:"
echo "  setup      Set up analysis"
echo "  results    Process results and generate HTML"
echo "  nan-repair Replace \`nan' values in PPGenotypeProbs.txt file with NA"
echo ""
exit "$1"
}

# Note that we use `"$@"' to let each command-line parameter expand to a 
# separate word. The quotes around `$@' are essential!
# We need TEMP as the `eval set --' would nuke the return value of getopt.
TEMP=$(getopt -o ha:sfc:o:m: --long action:,submit,force,config:,output:,main-directory:,help,data: \
     -n 'example.bash' -- "$@")

if [[ $? != 0 ]]; then echo "Terminating..." >&2 ; exit 1 ; fi

# Note the quotes around `$TEMP': they are essential!
eval set -- "$TEMP"

while true ; do
        case "$1" in
                --submit)
			echo "Submit"
			SUBMIT="yes"
			shift ;;
                -h|--help)
			usage 0
			shift ;;
                -f|--force)
			echo "Force"
			FORCE="yes"
			shift ;;
		-c|--config)
			ANALYSES_CONF="$2"
			if [[ ! -r "${ANALYSES_CONF}" ]]
			then
				echo "Can't read file '${ANALYSES_CONF}'."
				exit 1
			else
				# echo "Config: ${ANALYSES_CONF}"
				echo -n
			fi
			shift 2 ;;

		-o|--output)
			RESULTS_DIR="$2"
			if [[ ! -r "${RESULTS_DIR}" ]]
			then
				echo "Can't read file '${RESULTS_DIR}'."
				exit 1
			else
				# echo "Output: ${RESULTS_DIR}"
				echo -n
			fi
			shift 2 ;;

		-m|--main-directory)
			MAIN_DIRECTORY="$2"
			if [[ ! -d "${MAIN_DIRECTORY}" ]]
			then
				echo "Directory '${MAIN_DIRECTORY}' doesn't exist."
				exit 1
			else
				# echo "Output: ${MAIN_DIRECTORY}"
				echo -n
			fi
			shift 2 ;;

		-a|--action)
			ACTION="$2"
			shift 2 ;;

		-d|--data)
			SOURCE_DATA_DIR="$2"
			shift 2 ;;

		--) shift ; break ;;
                *) echo "Internal error! '$1'" ; exit 1 ;;
        esac
done
# echo "Remaining arguments:"
# for arg do echo '--> '"\`$arg'" ; done

if [[ -z "${ACTION}" ]]
then
	echo "No action has been given."
	usage 1
fi

if [[ -z "${SOURCE_DATA_DIR}" ]]
then
	echo "No data directory given."
	usage 1
fi

if [[ -z "${ANALYSES_CONF}" ]]
then
	echo "No config file has been given."
	usage 1
fi

function relative_to_absolute_dir() {
pushd "$1" > /dev/null
echo -n "$(pwd)"
popd > /dev/null
}

function action_results() {
	DIR="$1"
	OUTPUT_DIR="$(relative_to_absolute_dir "$2")"
 	pushd "../${DIR}/hapmap" > /dev/null
	# Convert ps to gif and copy the files
	for FILE in $(find . -name '*.ps')
	do
		# echo "Processing '${FILE}'"
		NEWFILE_BASE="$(echo -n "${FILE}" | sed -e 's/\.ps$//')"
		NEWFILE="${NEWFILE_BASE}.gif"
		NEWFILE_THUMB="${NEWFILE_BASE}-thumb.gif"
		if [[ ! -r "${NEWFILE}" ]]
		then
			echo "Converting '${FILE}' to '${NEWFILE}'"
			convert "${FILE}" -rotate 90 "${NEWFILE}"
			convert "${FILE}" -rotate 90 -resize 240x240 "${NEWFILE_THUMB}"
		fi
	done
	popd > /dev/null
}

while read LINE
do
	# Strip comments
	if [[ "$(echo ${LINE} | head -c 1)" = "#" ]]; then continue; fi
	DIR=$(echo ${LINE} | cut -d: -f1)
	OPTIONS=$(echo ${LINE} | cut -d: -f2)

	# 1. Copy the skeleton
	case "$ACTION" in
	setup)
		echo "1. Making the directory."
		# Create the directory (inside MAIN_DIRECTORY)
		if [[ "${PURGE_DIR}" = "yes" ]]
		then
			if [[ -r "../${MAIN_DIRECTORY}/${DIR}" && "${FORCE}" = "no" ]]
			then
				echo "../${MAIN_DIRECTORY}/${DIR} exists; please delete it by hand to proceed."
				exit 1
			elif [[ -r "../${MAIN_DIRECTORY}/${DIR}" && "${FORCE}" = "yes" ]]
			then
				rm -rf "../${MAIN_DIRECTORY}/${DIR}"
				rm -rf "../${DIR}"
			fi
		fi
		# Directory for results
		mkdir -p -v "${HOME}/${MAIN_DIRECTORY}/${DIR}"
		# Linking the directory to home
		ln -sv "${HOME}/${MAIN_DIRECTORY}/${DIR}" "../${DIR}"
		if [ "$?" != "0" ]
		then
			echo "Linking has failed."
			exit 1
		fi

		# Copy using tar to handle .bzr directory exclusion
		tar cf - . | tar xf - -C ../${DIR} --exclude .bzr --exclude html
		pushd "../${DIR}"

		# 2. Place the data
		echo "2. Placing the data."
		pushd hapmap
		for POPULATION in ${RESDIRS}
		do
			if [[ -r "../${POPULATION}" ]]
			then
				echo "${POPULATION} exists, please remove it by hand to proceed."
				exit 1
			fi
			mkdir -v -p "${POPULATION}"
			if [[ -L "${POPULATION}/chr22data" ]]
			then
				echo "${POPULATION}/chr22data exists already, removing."
				rm "${POPULATION}/chr22data"
			fi
			ln -s "${SOURCE_DATA_DIR}/${POPULATION}/chr22data" "${POPULATION}/chr22data"
			RES="$?"
			if [[ "${RES}" != "0" ]]
			then
				echo "ln failed: ${RES}"
				tree
				exit ${RES}
			fi
		done
		popd

		# 3. Set the directory in the PBS script
		echo "3. Setting up pbs script"
		echo "PBS_SCRIPT(${HOME}/${DIR}, ${DIR})" | m4 pbs.m4 - > "hapmap/${PBS_FILE}"

		# 4. Write the options (how?!)
		echo "4. Writing options to file."
		echo "${OPTIONS}" > hapmap/extra-cli-options.txt

		# 5. Submit the task
		# echo qsub "${PBS_FILE}"
		if [[ "${SUBMIT}" = "yes" ]]
		then
			echo "5. Submitting the task to PBS."
			pushd hapmap
			qsub "${PBS_FILE}"
			popd
		else
			echo "Not submitting cluster"
		fi
		popd
		;;
	results)
		# rm -rf "${RESULTS_DIR}"
		# mkdir -p "${RESULTS_DIR}"
		action_results "${DIR}" "${RESULTS_DIR}"

		# Results are now being processed separately
		;;
	nan-repair)
		cp -v hapmap/MutualInformation.R "../${DIR}/hapmap"
		cp -v hapmap/calculate-mi.sh "../${DIR}/hapmap"
		pushd "../${DIR}/hapmap"
		for POPULATION in ${RESDIRS}
		do
			echo "nan(${POPULATION})"
			sed -i -e 's/nan/NA/g' "${POPULATION}/${HAPMIXMAP_RES_SUBDIR}/PPGenotypeProbs.txt"
		done
		bash calculate-mi.sh
		popd
		;;
	*)
		echo "Unknown action ${ACTION}"
		exit 1
		;;
	esac
done < "${ANALYSES_CONF}"

if [[ "${ACTION}" = "results" ]]
then
	# write_index
	# python generate-html.py "${DIR}" "${OUTPUT_DIR}" "../${MAIN_DIRECTORY}"
	python generate-html.py --config "${ANALYSES_CONF}" --output "${OUTPUT_DIR}" --main-dir "${MAIN_DIRECTORY}"
fi

