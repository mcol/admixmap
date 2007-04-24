#!/usr/bin/perl
use strict; 
use File::Path;

my $function_file = "../../doanalysis.pl";

require $function_file or die("cannot find doanalysis.pl");

# Change this to the location of the admixmap executable
my $executable = '../../admixmap';

# Change this to the location of the R script
my $rscript = "../../AdmixmapOutput.R";

# $arg_hash is a hash of parameters passed to
# the executable as arguments.
#
# keys (left-hand side) are parameter names
# values (right-hand side) are parameter values
my $arg_hash = 
{
    samples                    => 2100, 
    burnin                     => 100,
    every                      => 1,
    displaylevel              => 0,
    locusfile                  => 'lociDRD2.txt',
    genotypesfile              => 'genotypesDRD2NAm.txt',
    populations                => 1,
    
# output files
    logfile                    => 'logfile.txt',
    hwscoretestfile            => 'HWtest.txt',
    resultsdir                 => 'DRD2freqresultsNAm',
    allelefreqoutputfile       => 'hapfreq.txt'
};

doAnalysis($executable, $rscript, $arg_hash);

$arg_hash->{resultsdir}             = 'DRD2freqresultsEur';
$arg_hash->{genotypesfile}             = 'genotypesDRD2Eur.txt';
doAnalysis($executable, $rscript, $arg_hash);

$arg_hash->{resultsdir}             = 'DRD2freqresultsAfr';
$arg_hash->{genotypesfile}             = 'genotypesDRD2Afr.txt';
doAnalysis($executable, $rscript, $arg_hash);

$arg_hash->{resultsdir}             = 'CandidatesfreqresultsEur';
$arg_hash->{locusfile}             = 'lociCandidates.txt';
$arg_hash->{genotypesfile}             = 'genotypesCandidatesEur.txt';
doAnalysis($executable, $rscript, $arg_hash);

$arg_hash->{resultsdir}             = 'CandidatesfreqresultsNAm';
$arg_hash->{genotypesfile}             = 'genotypesCandidatesNAm.txt';
doAnalysis($executable, $rscript, $arg_hash);



