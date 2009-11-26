#!/bin/bash


# Rename loci in the locus file from "999" to "X999", which matches the names
# used in the genotype files:
sed -i 's/^[0-9]/X\0/' data/loci.txt

# Re-sort the prior allele frequencies:
cut -c 2- data/priorallelefreqs.txt \
	| sort -n \
	| sed 's/^/X/;s/^Xoc/loc/' > data/x \
    && mv data/x data/priorallelefreqs.txt
