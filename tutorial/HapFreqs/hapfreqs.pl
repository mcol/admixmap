#!/usr/bin/perl
use strict; 
use File::Path;

# Change this to the location of the admixmap executable
my $executable = '../../test/admixmap.exe';

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
    coutindicator              => 0,
    locusfile                  => 'lociDRD2.txt',
    genotypesfile              => 'genotypesDRD2NAm.txt',
    populations                => 1,
	analysistypeindicator  => 1, 
    
# output files
    logfile                    => 'logfile.txt',
    hwscoretestfile            => 'HWtest.txt',
    resultsdir                 => 'DRD2freqresultsNAm',
    allelefreqoutputfile       => 'hapfreq.txt'
    };

doAnalysis($executable,$arg_hash);

$arg_hash->{resultsdir}             = 'DRD2freqresultsEur';
$arg_hash->{genotypesfile}             = 'genotypesDRD2Eur.txt';
doAnalysis($executable,$arg_hash);

$arg_hash->{resultsdir}             = 'DRD2freqresultsAfr';
$arg_hash->{genotypesfile}             = 'genotypesDRD2Afr.txt';
doAnalysis($executable,$arg_hash);

$arg_hash->{resultsdir}             = 'CandidatesfreqresultsEur';
$arg_hash->{locusfile}             = 'lociCandidates.txt';
$arg_hash->{genotypesfile}             = 'genotypesCandidatesEur.txt';
doAnalysis($executable,$arg_hash);

$arg_hash->{resultsdir}             = 'CandidatesfreqresultsNAm';
$arg_hash->{genotypesfile}             = 'genotypesCandidatesNAm.txt';
doAnalysis($executable,$arg_hash);


sub getArguments
{
    my $hash = $_[0];
    my $arg = '';
    foreach my $key (keys %$hash){
	$arg .= ' --'. $key .'='. $hash->{$key};
    }
    return $arg;
}

sub doAnalysis
{
    my ($prog,$args) = @_;
    my $command = $prog.getArguments($args);
    if (-e $args->{resultsdir}) {
	rmtree($args->{resultsdir});
    } 
    mkpath($args->{resultsdir});
    $ENV{'RESULTSDIR'} = $args->{resultsdir};
    print "Results will be written to subdirectory $ENV{'RESULTSDIR'}\n";
    system($command);
    print "Starting R script to process output\n";
    system("RCMD BATCH --quiet --no-save --no-restore ..\\..\\test\\AdmixmapOutput.R $args->{resultsdir}/Rlog.txt");
    print "R script completed\n\n";
}

