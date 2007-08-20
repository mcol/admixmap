#!/bin/bash

POP=ceu

for i in $(seq 1 22)
do
	ln -s genotypes-chr${i}-${POP}-r21-nr-fwd-legend.txt chr${i}_legend.txt
	ln -s genotypes-chr${i}-${POP}-r21-nr-fwd-phased chr${i}_phased.txt
	ln -s genotypes-chr${i}-${POP}-r21-nr-fwd-sample.txt chr${i}_sample.txt
done
