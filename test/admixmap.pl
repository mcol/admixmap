#!/usr/bin/perl -w
use strict;

my $DEBUG = 0; # zero gives less output

# Change this to the location of the admixmap executable
my $executable = './admixmap';
my $resultsdir = "results";

# $arg_hash is a hash of parameters passed to
# the executable as arguments.
#
# keys (left-hand side) are parameter names
# values (right-hand side) are parameter values
my $arg_hash = 
{
    stratificationtestfile       => 'strat_test.txt',
    logfile                      => 'log_test.txt',
    genotypesfile                => 'data/genotypes.txt',
    locusfile                    => 'data/loci.txt',
    priorallelefreqfile          => 'data/priorallelefreqs.txt',
    outcomevarfile               => 'data/outcomevar_diabetes.txt',
    samples  => 10,#22000,
    burnin   => 5,#2000,
    every    => 1,

resultsdir => "$resultsdir",
    paramfile                    => 'paramfile.txt',
    regparamfile => 'regression.txt',
    indadmixturefile             => 'ind_admixture.txt',
    ergodicaveragefile           => 'ergodicaverage.txt',

    allelicassociationscorefile  => 'locusscore_test.dat',
    ancestryassociationscorefile => 'locusscore_test2.dat',
    affectedsonlyscorefile        => 'affectedsonlyscorefile.txt',
    analysistypeindicator     => 3,
    coutindicator   => 1,
    targetindicator => 0,
    covariatesfile               => 'data/covariates3.txt',
    haplotypeassociationscorefile => 'hapassocscore.txt',
    allelefreqoutputfile         => 'allelefreqoutputfile_test.dat'

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
    print "Starting R script to process output\n";
    system('RCMD BATCH --quiet --no-save --no-restore AdmixmapOutput.R results/Rlog.txt');
    print "R script completed\n\n";
}



