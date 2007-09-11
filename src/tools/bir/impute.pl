#!/usr/bin/perl
use strict;

my $dataprefix = "data/chr22_5kloci/impute";

my @Panels = ("CEU", "YRI", "JPTCHB");
my @Ne      = (11418, 17469, 14269);#from the manual

#outdp = number of decimal places in output
#os 1 means output loci in haplotype file but not genotype file, ie the masked loci
my $options =  "-outdp 6 -os 1";


for (my $i=0; $i < 3; ++$i){
    my $datadir = "$dataprefix/@Panels[$i]";
    my $resultsdir = "results/impute/@Panels[$i]";

    system( "impute -h $datadir/haplo.txt -l $datadir/legend.txt -g $datadir/geno.txt -m $datadir/map.txt -Ne @Ne[$i] -i $resultsdir/info.txt -o $resultsdir/out.txt $options\n\n");
};

