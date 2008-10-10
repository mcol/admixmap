#!/usr/bin/perl
use strict;
use File::Path;
use Cwd;

my $dir = getcwd;

my @Panels = ("CEU", "YRI", "JPTCHB");
my $states = 8; 

# -T number of random starts
# -C number of EM iterations
# -K number of states
# -S random seed (default is from system time)
# -s sample phased haplotypes 
# -p print parameter estimates
# -M2 fixed mixture proportions
# -m print estimated probabilities for missing genotypes - causes crash  
# -H-4 turn off haplotype inference - causes empty output file
# -H<n> for n samples per start 
# -o output prefix
# -b haplotype file
#my $options =  "-T2 -C50 -K$states -M2 -m"; 
my $options =  "-T20 -C50 -K$states -M2 -s10000 -H10 -S7001"; 

my $whichchr = "chr22_5kloci";
#my $whichchr = "chr22";

my $dataprefix = "data/$whichchr/fastphase";
my $mytempdir = "/exports/work/scratch/pmckeigu";

for my $pop (@Panels) {
    my $datadir = "$dataprefix/$pop";
    my $resultsdir = "$mytempdir/$whichchr/fastphase/$pop";
    if(!(-e "$resultsdir")) {
 	 mkpath "$resultsdir" or die "cannot make directory $resultsdir";
 	 }
    my $command = "fastPHASE $options -b$datadir/fphaplotypes.inp -o'$TMPDIR'/fastphase$pop $datadir/fastphase.inp";
    open(SCRIPT, "> $datadir/fastphase.sh"); # shell script in 
    print(SCRIPT "cd $dir\n"); # set working directory
    print(SCRIPT " $command\n"); 
	close(SCRIPT);
    system("qsub -cwd -e $dir/fastphase$poperr.txt -o $dir/fastphase$popout.txt -l h_rt=01:00:00 -N fastphase$pop $datadir/fastphase.sh");
};
