#!/usr/bin/perl -w
use strict; 
use File::Path;

my $executable = './admixmap';

my $arg_hash = {
    burnin   => 10, 
    samples  => 60,
    every    => 1,
    locusfile                    => "IndData/loci.txt",
    genotypesfile                => "IndData/genotypes.txt",
    priorallelefreqfile          => "IndData/priorallelefreqs3way.txt",
    fixedallelefreqs => 1,
    randommatingmodel            => 1,
    globalrho                    => 0,
    initalpha0                   => "1,1,1", # parameter vectors for Dirichlet prior on admixture 
    initalpha1                   => "1,1,0",
    chib                         => 1,
    thermo                       => 1,
    numannealedruns             => 100,

    resultsdir => "IndResults",
    logfile                      => "logfile.txt",
    indadmixturefile             => "indadmixture.txt"
};

## no admixture model: one Afr, one Eur parent
$arg_hash->{initalpha0} = "0,1,0";
$arg_hash->{initalpha1} = "1,0,0";
$arg_hash->{resultsdir} = "IndResults010-100";
#&doAnalysis($executable, $arg_hash);

## simplest admixture model: one Afr/Eur, one Afr parent 
$arg_hash->{initalpha0} = "1,0,0"; 
$arg_hash->{initalpha1} = "1,1,0";
$arg_hash->{resultsdir} = "IndResults010-110";
&doAnalysis($executable, $arg_hash);

sub getArguments {
    my $hash = $_[0];
    my $arg = '';
    foreach my $key (keys %$hash){
	$arg .= ' --'. $key .'='. $hash->{$key};
    }
    return $arg;
}

sub doAnalysis {
    my ($prog,$args) = @_;
    my $command = $prog.getArguments($args);
    if (-e $args->{resultsdir}) {
	rmtree($args->{resultsdir});
    } 
    mkpath($args->{resultsdir});
    $ENV{'RESULTSDIR'} = $args->{resultsdir};
    print "\nResults will be written to subdirectory $ENV{'RESULTSDIR'}\n";
    system($command);
    my $rcmd = "R CMD";
    if($^O eq "MSWin32") {
	$rcmd = "Rcmd";
    }
    print "Starting R script to process output\n";
    system("$rcmd BATCH --quiet --no-save --no-restore ../test/AdmixmapOutput.R $args->{resultsdir}/Rlog.txt\n");
    print "R script completed\n\n";
}

