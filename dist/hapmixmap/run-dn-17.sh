#!/bin/bash
# 
# This runs hapmixmap with tutorial configuration files, with the
# diabetic nephropathy data set, containing 17 typed loci.

hapmixmap training-initial.conf
hapmixmap training-resume.conf
hapmixmap -ftesting.conf \
	--ccgenotypesfile=data/CaseControlGenotypes_2_genotypes.txt \
	--outcomevarfile=data/CaseControlGenotypes_2_outcome.txt
