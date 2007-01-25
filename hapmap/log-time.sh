#!/bin/bash

# Logs time usage in format:
#
# Population, States, Samples, Haploid, Diploid, Elapsed
#
# Given argument: File with hapmixmap output.

LOG_FILE="$1"
TIME_FILE=time-stats.txt

POPULATION=$(cat $LOG_FILE | grep Loading | head -n 1 | cut -d" " -f2 | cut -d"/" -f1)
STATES=$(cat $LOG_FILE | grep states | cut -d" " -f 5)
SAMPLES=$(cat $LOG_FILE | grep iterations | cut -d" " -f1)
HAPLOID=$(cat $LOG_FILE | grep haploid | cut -d" " -f4)
DIPLOID=$(cat $LOG_FILE | grep diploid | cut -d" " -f1)
ELAPSED=$(cat $LOG_FILE | grep "Elapsed seconds" | cut -d"=" -f2)

echo -e >> "$TIME_FILE" "$POPULATION\t$STATES\t$SAMPLES\t$HAPLOID\t$DIPLOID\t$ELAPSED"
