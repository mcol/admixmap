#!/usr/bin/perl
use strict;
use File::Path;

my $whichchr = "chr22_5kloci";

my$dataprefix = "data/$whichchr/fastphase";

my @Panels = ("YRI", "CEU", "JPTCHB");
my $states = 6;
my $options =  "-T2 -C50 -K$states -s1000 -p -M2";
# -T=number of random starts
# -C=number of EM iterations
# -K=number of states
# -S=random seed (default is from system time)
# -p=print parameter estimates
#-M2=fixed mixture proportions
#-m=print estimated probabilities for missing genotypes
#-o=output files prefix
#-b=haplotype file

for (my $i = 0;  $i < 3; ++$i) {
    my $datadir = "$dataprefix/@Panels[$i]";
    my $resultsdir = "results/$whichchr/fastphase/@Panels[$i]/";
    if(!(-e "$resultsdir")) { 
	mkpath "$resultsdir" or die "cannot make directory $resultsdir";
    }
    system("fastPHASE $options -b$datadir/haplotypes.inp -o$resultsdir $datadir/fastphase.inp");
};

