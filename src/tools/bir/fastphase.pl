#!/usr/bin/perl
use strict;
use File::Path;

my $dataprefix = "data/chr22_5kloci/fastphase";

my @Panels = ("CEU", "YRI", "JPTCHB");
my $states = 6;

# -T=number of random starts
# -C=number of EM iterations
# -K=number of states
# -S=random seed (default is from system time)
# -p=print parameter estimates
#-M2=fixed mixture proportions
#-m=print estimated probabilities for missing genotypes
#-o=output prefix
#-b=haplotype file
my $options =  "-T20 -C50 -K$states -S1000 -p -M2 ";


for my $pop (@Panels) {
    my $datadir = "$dataprefix/$pop";
    my $resultsdir = "results/fastphase/$pop";
    if(!(-e "$resultsdir")) {
 	 mkpath "$resultsdir" or die "cannot make directory $resultsdir";
 	 }	
    system("fastPHASE $options -b$datadir/fphaplotypes.inp -o$resultsdir/fastphase $datadir/fastphase.inp");
};

