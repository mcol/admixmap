#!/usr/bin/perl
use strict; 
use File::Path;

################### DO NOT EDIT ABOVE THIS LINE ########################
# Change this to the name of the program to be run
my $executable = 'c:/cvs/genepi/test/admixmap';

my $arg_hash = 
{
    burnin   => 50,
    samples  => 250,
    every    => 1,
    #numannealedruns => 0,
    indadmixhiermodel            => 0,
    #chib                         => 1,
    randommatingmodel            => 1,
    globalrho                    => 0,
    fixedallelefreqs             => 1,
    locusfile                    => 'data/loci.txt',
    genotypesfile                => 'data/genotypes.txt',
    priorallelefreqfile          => 'data/trueallelefreqs.txt',
    logfile                      => 'log.txt',
    indadmixturefile             => 'individualVarSamples.txt',
    #indadmixmodefile             => 'posteriormode.txt'
};

$arg_hash->{resultsdir}         = 'results-11-11';
$arg_hash->{admixtureprior}         = "1,1";
$arg_hash->{admixtureprior1}        = "1,1";
doAnalysis($executable,$arg_hash);

$arg_hash->{resultsdir}         = 'results-10-11';
$arg_hash->{admixtureprior}         = "0,1";
$arg_hash->{admixtureprior1}        = "1,1";
#doAnalysis($executable,$arg_hash);

sub getArguments {
    my $hash = $_[0];
    my $arg = '';
    foreach my $key (keys %$hash) {
	$arg .= ' --'. $key .'='. $hash->{$key};
    }
    return $arg;
}

sub doAnalysis {
    my ($prog, $args) = @_;
    my $command = $prog.getArguments($args);
    if (-e $args->{resultsdir}) {
	rmtree($args->{resultsdir});
    } 
    mkpath($args->{resultsdir});
    $ENV{'RESULTSDIR'} = $args->{resultsdir};
    print "Results will be written to subdirectory $ENV{'RESULTSDIR'}";
    system($command);
    print "Starting R script to process output\n";
    system("R CMD BATCH --quiet --no-save --no-restore \
           c:/cvs/genepi/test/AdmixmapOutput.R $args->{resultsdir}/Rlog.txt");
    print "R script completed\n\n";
}

