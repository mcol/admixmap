#!/usr/bin/perl
use strict;
use File::Path;

my $whichchr = "chr22_5kloci";

my$dataprefix = "data/$whichchr/fastphase";

my @Panels = ("YRI", "CEU", "JPTCHB");
my $states = 8;
#my $options =  "-T2 -C50 -K$states -M2 -m";
my $options =  "-T10 -C50 -K$states -M2 -s1000";
# -T=number of random starts
# -C=number of EM iterations
# -K=number of states
# -S=random seed (default is from system time)
# -s=number of samples of phased haplotypes
# -p=print parameter estimates
# -M2=fixed mixture proportions

# -m=print estimated probabilities for missing genotypes - causes crash 
#? not in manual for v1.2

# -H-4 = turn off haplotype inference - causes empty output file
# -o=output files prefix
# -b=haplotype file

for (my $i = 0;  $i < 3; ++$i) {
    my $datadir = "$dataprefix/@Panels[$i]";
    my $resultsdir = "results/$whichchr/fastphase/@Panels[$i]";
    if(!(-e "$resultsdir")) { 
	mkpath "$resultsdir" or die "cannot make directory $resultsdir";
    }
    system("fastPHASE $options -b$datadir/fphaplotypes.inp -o$resultsdir/ $datadir/fastphase.inp");
};

