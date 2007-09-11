#!/usr/bin/perl
use strict;
use File::Path;

my $whichchr = "chr22_5kloci";
#my $whichchr = "chr22";

my$dataprefix = "data/$whichchr/impute";

my @Panels = ("CEU", "YRI", "JPTCHB");
my @Ne      = (11418, 17469, 14269);#from the manual

#outdp = number of decimal places in output
#os 1 means output loci in haplotype file but not genotype file, ie the masked loci
my $options =  "-outdp 6 -os 1";


for (my $i=0; $i < 3; ++$i){
    my $datadir = "$dataprefix/@Panels[$i]";
    my $resultsdir = "results/impute/@Panels[$i]";
    if(!(-e "$resultsdir")) {
 	 mkpath "$resultsdir" or die "cannot make directory $resultsdir";
    }

    system( "impute -h $datadir/haplo.txt -l $datadir/legend.txt -g $datadir/geno.txt -m $datadir/map.txt -Ne @Ne[$i] -i $resultsdir/info.txt -o $resultsdir/out.txt $options\n\n");
};

