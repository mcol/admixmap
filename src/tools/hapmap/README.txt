This directory contains:

* FUD, a program to format unphased (diploid) HapMap data for HAPMIXMAP

Formatting unphased HapMap data
---------------------------------
Note: This program is obsolete and unsupported. HAPMIXMAP is intended to be used with phased HapMap data. This section is for information only.

You will have to download the unphased data yourself. For each population, you will need the information file, pedinfo2sample_XXX.txt.gz (where XXX is a the population code), currently located at http://www.hapmap.org/downloads/samples_individuals/.
Then, for each chromosome, you will need the 'non-redundant' fwd strand genotypes file, genotypes_chrY_XXX_r22_nr.bZZ_fwd.txt.gz (where Y is the chromosome number and ZZ is the current build).
You will then need to run the perl script repairSecondColumn.pl to ensure the foramatting program reads the bases correctly. Finally run the formatting program, FUD. Type FUD -h for a list of options.




