#!/usr/bin/perl
use strict;
use File::Path;

my$dataprefix = "data/chr22_5kloci/fastphase";

my @Panels = ("CEU", "YRI", "JPTCHB");
my $states = 6;
my $options =  "-T20 -C50 -K$states -s1000 -p -M2";


for my $pop in @Panels {
    my $datadir = "$dataprefix/$pop";
    my $resultsdir = "results/$whichchr/fastphase/@Panels[$i]/";
    if(!(-e "$resultsdir")) { 
	mkpath "$resultsdir" or die "cannot make directory $resultsdir";
    }
    # -b haplotypes file
    # -o prefix to results file names
    system("fastphase $options -b$datadir/haplotypes.inp -o$resultsdir $datadir/fastphase.inp");
};

