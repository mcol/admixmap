#!/usr/bin/perl
use strict;
use File::Path;
use Cwd;

my $dir = getcwd;

my @Panels = ("CEU", "YRI", "JPTCHB");
my $states = 8; 

# -T=number of random starts
# -C=number of EM iterations
# -K=number of states
# -S=random seed (default is from system time)
# -s=number of samples of phased haplotypes 
# -p=print parameter estimates
#-M2=fixed mixture proportions
# -m=print estimated probabilities for missing genotypes - causes crash  
# -H-4 = turn off haplotype inference - causes empty output file 
#-o=output prefix
#-b=haplotype file
#my $options =  "-T2 -C50 -K$states -M2 -m"; 
my $options =  "-T20 -C25 -K$states -M2 -s1000"; 

my $whichchr = "chr22_5kloci";
#my $whichchr = "chr22";

my$dataprefix = "data/$whichchr/fastphase";

for my $pop (@Panels) {
    my $datadir = "$dataprefix/$pop";
    my $resultsdir = "results/$whichchr/fastphase/$pop";
    if(!(-e "$resultsdir")) {
 	 mkpath "$resultsdir" or die "cannot make directory $resultsdir";
 	 }
    my $command = "fastPHASE $options -b$datadir/fphaplotypes.inp -o$resultsdir/fastphase $datadir/fastphase.inp";
    open(SCRIPT, "> $datadir/fastphase.sh");
    print(SCRIPT "cd $dir\n");
    print(SCRIPT " $command\n");
	close(SCRIPT);
    system("qsub -l h_rt=01:00:00 $datadir/fastphase.sh");
};


