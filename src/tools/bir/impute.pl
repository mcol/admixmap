#!/usr/bin/perl
use strict;

my$dataprefix = "data/chr22_5kloci/impute";

my @Panels = ("YRI", "CEU", "JPTCHB");
my $options =  "-Ne 11400 -outdp 6";

for my $pop(@Panels) {
    my $datadir = "$dataprefix/$pop";
    system("impute -h $datadir/haplo.txt -l $datadir/legend.txt -g $datadir/geno.txt -m $datadir/map.txt -o impgenotypes$pop.txt $options");
};

