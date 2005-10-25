#!/usr/bin/perl -w
use strict;
use File::Path;

# Change this to the location of the admixmap executable
my $executable = '../test/admixmap';

# $arg_hash is a hash of parameters passed to
# the executable as arguments.
#
# keys (left-hand side) are parameter names
# values (right-hand side) are parameter values
my $arg_hash = 
{
    samples                    => 1200, 
    burnin                     => 200,
    every                      => 5,
    locusfile                  => 'data/x_loci.txt',
    genotypesfile              => 'data/genotypes.txt',
    xonlyanalysis              => 1,
    coutindicator              => 1,
    priorallelefreqfile        => 'data/priorallelefreqs_x.txt',
    resultsdir                 => 'test',  
    
# output files
    logfile                    => 'logfile.txt',
    paramfile                  => 'param.txt',
    ergodicaveragefile         => 'ergodicaverage.txt',
     
# extra output files
    indadmixturefile           => 'indadmixture.txt'
    };

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
    system("RCMD BATCH --quiet --no-save --no-restore ../test/AdmixmapOutput.R $args->{resultsdir}/Rlog.txt");
    print "R script completed\n\n";
}

