#!/bin/bash

# Generates configuration file for analyses.

if [ -z "$1" ]
then
	ANALYSES_CONF="analyses.txt"
else
	ANALYSES_CONF="$1"
fi

if [ -r "${ANALYSES_CONF}" ]
then
	echo "File '${ANALYSES_CONF}' is present. Please delete it by hand to proceed."
	exit 1
fi

cat lib-analyses.txt | grep -v -E "^#" | tee "${ANALYSES_CONF}"

# http://actin.ucd.ie/trac/genepi/ticket/85
for SEED in $(seq 1 4); do
	# Arrival rate prior
	for ARP in \
		"1.2, 0.05, 1, 0.5" \
		"12,  0.5,  1, 0.5" \
		"1.2, 0.05, 2" \
		"12,  0.5,  2"
	do
	# for ARP in "12, 0.5, 2"; do
		# Allele frequency precision prior
		for AFPP in "2, 8" "0.25, 1"; do
			for MIXPROP in 1
			do
				# DIR_NAME="coc-4-fmp-${MIXPROP}-s${SEED}"
				DIR_NAME="irw-2-arp-${ARP}-afpp-${AFPP}-s${SEED}"

				# Slugify the directory name
				DIR_NAME="$(echo -n "${DIR_NAME}" | sed -e 's/[^[:alnum:]\.]/-/g' | tr -s '-')"
				echo -n "${DIR_NAME}:"
				# echo -n "--mixturepropsprecisionprior=\"${MIXPPP}\" "
				# echo -n "--fixedmixtureprops=0 --fixedmixturepropsprecision=0 "
				echo -n " --arrivalrateprior=\"${ARP}\""
				echo -n " --allelefreqprecisionprior=\"${AFPP}\""
				echo -n " --fixedmixtureprops=\"${MIXPROP}\""
				if [[ "${MIXPROP}" == "0" ]]
				then
					echo -n " --mixturepropsprecision=\"20\""
				fi
				echo -n " --seed=\"${SEED}\""
				echo
			done | tee -a "${ANALYSES_CONF}"
		done
	done
done
