#!/usr/bin/perl
use strict;

my$dataprefix = "data/chr22_5kloci/fastphase";

my @Panels = ("CEU", "YRI", "JPTCHB");
my $states = 6;
my $options =  "-T20 -C50 -K$states -s1000 -p -M2";


for my $pop in @Panels {
    my $datadir = "$dataprefix/$pop";
    system("fastphase $options -bfphaplotypes.inp fastphase.inp ");
};

