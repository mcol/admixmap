#!/usr/bin/perl -w
use strict;

my $DEBUG = 0; # zero gives less output

# Change this to the location of the admixmap executable
my $executable = './admixmap';

# $arg_hash is a hash of parameters passed to
# the executable as arguments.
#
# keys (left-hand side) are parameter names
# values (right-hand side) are parameter values
my $arg_hash = 
{
#data files
    genotypesfile                   => 'sim/genotypes.txt',
    locusfile                          => 'sim/loci.txt',
    priorallelefreqfile             => 'sim/priorallelefreqs.txt',
    # covariatesfile                  => 'data/covariates3.txt',
    outcomevarfile               => 'sim/outcome.txt',

#main options
    analysistypeindicator     => 2, # continuous outcome 
    coutindicator   => 1, #verbose output
    samples  => 500,
    burnin   => 100,
    every    => 5,

#output files
    resultsdir               => "sim",
    logfile                     => 'logfile.txt',
    paramfile               => 'paramfile.txt',
    indadmixturefile     => 'indadmixture.txt',
    ergodicaveragefile => 'ergodicaverage.txt'
    # allelefreqoutputfile  => 'allelefreqoutputfile.dat',

#optional tests
    # allelicassociationscorefile       => 'allelicassociationscorefile.dat',
    # ancestryassociationscorefile  => 'ancestryassociationscorefile.dat',
    # affectedsonlyscorefile             => 'affectedsonlyscorefile.txt',
    # haplotypeassociationscorefile => 'hapassocscore.txt',
    # stratificationtestfile                   => 'strat_test.txt'
};

#doAnalysis($executable,$arg_hash);

$arg_hash->{outcomevarfile} = 'sim/outcome.txt';
$arg_hash->{analysistypeindicator} = 3;
$arg_hash->{regparamfile} = 'regparamfile.txt';
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
    unless (-e "results"){
	system("mkdir results");
    }
    print $command if $DEBUG;
    system($command);
    
# Comment out the next three lines to run admixmap without R script
    print "Starting R script to process output\n";
    system("R --quiet --no-save --no-restore <AdmixmapOutput.R >results/Rlog.txt RESULTSDIR=$arg_hash->{resultsdir}");
    print "R script completed\n\n";
}



