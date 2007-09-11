#!/usr/bin/perl
use strict;

my$dataprefix = "data/chr22_5kloci/impute";

my @Panels = ("CEU", "YRI", "JPTCHB");
my @Ne      = (11418, 17469, 14269);#from the manual

#outdp = number of decimal places in output
#os 2 means output only loci in both genotype file and haplotype file, ie the masked loci
my $options =  "-outdp 6 -os 2";


for (my $i=0; $i < 3; ++i){
    my $datadir = "$dataprefix/@Panels[$i]";
    my $resultsdir = "results/impute/@Panels[$i];

    system("impute -h $datadir/haplo.txt -l $datadir/legend.txt -g $datadir/geno.txt -m $datadir/map.txt -o impgenotypes$pop.txt -Ne@Ne[$i] -i $resultsdir/info.txt -o $resultsdir/out.txt $options");
};

