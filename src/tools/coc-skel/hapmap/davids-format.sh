#!/bin/bash

for i in $(seq 1 22)
do
	# -c=foo <- yuck!
	FPHD -c=${i} -p=hapmap-source -g=pddata/gt${i}.txt -l=pddata/loci${i}.txt -v
done

