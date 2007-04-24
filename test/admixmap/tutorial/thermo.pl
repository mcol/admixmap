#!/usr/bin/perl
use strict; 
use File::Path;
use Parallel::MPI::Simple;

my $function_file = "../doanalysis.pl";

require $function_file or die("cannot find doanalysis.pl");

# Change this to the location of the admixmap executable
my $executable = '../admixmap';

# Change this to the location of the R script
my $rscript = "../AdmixmapOutput.R";

# command-line options are stored in an associative array (known as a hash in perl)  
my $arg_hash = {
#data files
    genotypesfile                   => 'data/genotypesAIMsOnly.txt',
    locusfile                       => 'data/lociAIMsOnly.txt',
#main options
    samples  => 25,
    burnin   => 5,
    every    => 1,
    thermo   => 1,
    numannealedruns => 50,  
    displaylevel => 2,
#output file options
    logfile                     => 'log.txt',
};


MPI_Init();
my $rank = MPI_Comm_rank(MPI_COMM_WORLD);
my $np = MPI_Comm_size(MPI_COMM_WORLD);
my $index = 0;

# model with reference prior on allele freqs in 1 population, skin reflectance as continuous outcome var
$arg_hash->{populations}           = 1;
$arg_hash->{resultsdir}            = 'SinglePopResults';
if($rank == $index % $np) {
    &doAnalysis($executable, $rscript, $arg_hash);
}
$index = $index + 1;

# model with reference prior on allele freqs in 2 populations
$arg_hash->{populations}           = 2;
#$arg_hash->{samples}   = 6000;
#$arg_hash->{burnin}    = 1000;
$arg_hash->{paramfile}                 = 'popadmixparams.txt',
$arg_hash->{resultsdir}            = 'TwoPopsResults';  
if($rank == $index % $np) {
    &doAnalysis($executable, $rscript, $arg_hash);
}
$index = $index + 1;

# model with reference prior on allele freqs in 3 populations
$arg_hash->{populations}           = 3;
$arg_hash->{resultsdir}            = 'ThreePopsResults';  
if($rank == $index % $np) {
    &doAnalysis($executable, $rscript, $arg_hash);
}



