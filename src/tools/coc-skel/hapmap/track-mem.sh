#!/bin/bash

if [ -z "$1" ]
then
	echo "Usage $0 procname outfile"
	exit 1
fi

if [ -z "$2" ]
then
	echo "Usage $0 procname outfile"
	exit 1
fi

PROCNAME="$1"

while true
do
	ps axv | grep "$PROCNAME" \
		| grep -v bash \
		| grep -v grep \
		| grep -v perl \
		| awk '{ print $7; }'
	sleep 1
done | tee "$2"
