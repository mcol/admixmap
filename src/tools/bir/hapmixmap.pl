#!/usr/bin/perl
use strict;
use File::Path;
use Cwd;

my $dir = getcwd;
my @Panels = ("CEU", "YRI", "JPTCHB");

my $whichchr = "chr22_5kloci";
#my $whichchr = "chr22";

my$dataprefix = "data/$whichchr/hapmixmap";

for my $pop (@Panels) {
    my $datadir = "$dataprefix/$pop";
    my $resultsdir = "results/$whichchr/hapmixmap/$pop";
    if(!(-e "$resultsdir")) {
 	 mkpath "$resultsdir" or die "cannot make directory $resultsdir";
 	 }
    my $command = "hapmixmap $configfile";
    open(SCRIPT, "> $datadir/hapmixmap.sh");
    print(SCRIPT "cd $dir\n");
    print(SCRIPT " $command\n");
	close(SCRIPT);
    system("qsub -l h_rt=01:00:00 $datadir/hapmixmap.sh");
};

