#!/usr/bin/perl -w
use strict;

my $DEBUG = 0; # zero gives less output

# Change this to the location of the admixmap executable
my $executable = './admixmap';
# Change this to desired name of results directory
my $resultsdir = "results";

# $arg_hash is a hash of parameters passed to
# the executable as arguments.
#
# keys (left-hand side) are parameter names
# values (right-hand side) are parameter values
my $arg_hash = 
{
#data files
    genotypesfile                   => 'data/genotypes.txt',
    locusfile                          => 'data/loci.txt',
    priorallelefreqfile             => 'data/priorallelefreqs.txt',
    covariatesfile                  => 'data/covariates3.txt',
    outcomevarfile               => 'data/outcomevar_diabetes.txt',

#main options
    analysistypeindicator     => 3,# binary outcome var
    coutindicator   => 0, #verbose output
    targetindicator => 0, # diabetes in column 1
    samples  => 20,
    burnin   => 5,
    every    => 1,

#output files
    resultsdir               => "$resultsdir",
    logfile                     => 'logfile.txt',
    paramfile               => 'paramfile.txt',
    regparamfile          => 'regparamfile.txt',
    indadmixturefile     => 'indadmixture.txt',
    ergodicaveragefile => 'ergodicaverage.txt',
    allelefreqoutputfile  => 'allelefreqoutputfile.dat',

#optional tests
    allelicassociationscorefile       => 'allelicassociationscorefile.dat',
    ancestryassociationscorefile  => 'ancestryassociationscorefile.dat',
    affectedsonlyscorefile             => 'affectedsonlyscorefile.txt',
    haplotypeassociationscorefile => 'hapassocscore.txt',
    stratificationtestfile                   => 'strat_test.txt'
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
    unless (-e "results"){
	system("mkdir results");
    }

    print $command if $DEBUG;
    system($command);

# Comment out the next three lines to run admixmap without R script
    print "Starting R script to process output\n";
    system("R --quiet --no-save --no-restore <AdmixmapOutput.R >results/Rlog.txt RESULTSDIR=$resultsdir");
    print "R script completed\n\n";
}



