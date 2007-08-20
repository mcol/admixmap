#!/bin/bash

# Logs time usage in format:
#
# Population, States, Samples, Haploid, Diploid, Elapsed
#
# Given argument: File with hapmixmap output.

LOG_FILE="$1"

function usage {
echo "Usage: $0 <log>"
exit $1
}

function debug_fields {
	for VAR in POPULATION STATES SAMPLES HAPLOID DIPLOID ELAPSED LOCI
	do
		eval "echo \"\$VAR: \${${VAR}}\""
	done
}

if [ ! -r "$LOG_FILE" ]
then
	echo "Can't read file '$LOGFILE'"
	usage 1
fi

TIME_FILE=time-stats.txt
CURRENT_TIME=$(date +%Y-%m-%d-%H-%M-%S)

POPULATION="$(cat $LOG_FILE | grep Loading | tail -n 1 | cut -d"/" -f2 | cut -d . -f 1)"
STATES="$(cat $LOG_FILE | grep states | cut -d" " -f 5)"
SAMPLES="$(cat $LOG_FILE | grep iterations | cut -d" " -f1)"
HAPLOID="$(cat $LOG_FILE | grep "haploid individuals" | cut -d" " -f4)"
if [ -z "$HAPLOID" ]
then
	HAPLOID=$(cat $LOG_FILE | grep "haploid individuals" | cut -d" " -f1)
fi
DIPLOID="$(cat ${LOG_FILE} | grep "diploid"| cut -d" " -f1)"
DIPLOID="${DIPLOID:=0}"
ELAPSED="$(cat ${LOG_FILE} | grep "Elapsed seconds" | cut -d"=" -f2)"
LOCI="$(cat "${LOG_FILE}" | grep "simple loci" | cut -d " " -f 1)"
BURNIN="$(cat "${LOG_FILE}" | grep "^burnin" | cut -d " " -f2 )"
AFTER_BURNIN=$(( $SAMPLES - $BURNIN ))

# echo -e >> "$TIME_FILE" "$CURRENT_TIME\t$POPULATION\t$STATES\t$SAMPLES\t$HAPLOID\t$DIPLOID\t$ELAPSED"
# cur_time	pop	states	samples	haplo	diplo	elaps	loci burnin after_burnin
echo -e "$CURRENT_TIME\t$POPULATION\t$STATES\t$SAMPLES\t$HAPLOID\t${DIPLOID}\t$ELAPSED\t${LOCI}\t${BURNIN}\t${AFTER_BURNIN}"
# echo -e "$POPULATION\t$STATES\t$SAMPLES\t$HAPLOID\t$DIPLOID\t$ELAPSED"
