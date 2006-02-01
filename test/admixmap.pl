#!/usr/bin/perl -w
use strict;

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
    genotypesfile                   => 'data/genotypes.txt',
    locusfile                          => 'data/loci.txt',
    priorallelefreqfile             => 'data/priorallelefreqs.txt',
    #covariatesfile                  => 'data/covariates2std.txt',
    #outcomevarfile               => 'data/outcomevars.txt',

#main options
    displaylevel   => 3, #verbose output
    #targetindicator => 0, # diabetes in column 1
    outcomes => 1,
    samples  => 600,
    burnin   => 100,
    every    => 5,
    numannealedruns  => 0,
    #indadmixhiermodel => 1,
    #fixedallelefreqs => 1,

#output files
    resultsdir               => "resultsNoAnneal",
    logfile                     => 'logfile.txt',
    paramfile               => 'paramfile.txt',
    # regparamfile          => 'regparamfile.txt',
    indadmixturefile     => 'indadmixture.txt',
    ergodicaveragefile => 'ergodicaverage.txt',
    allelefreqoutputfile  => 'allelefreqoutputfile.dat',

#optional tests
    #allelicassociationscorefile       => 'allelicassociationscorefile.dat',
    #ancestryassociationscorefile  => 'ancestryassociationscorefile.dat',
    #affectedsonlyscorefile             => 'affectedsonlyscorefile.txt',
    #haplotypeassociationscorefile => 'hapassocscore.txt',
    stratificationtestfile                   => 'strat_test.txt'
};


#doAnalysis($executable,$arg_hash);

$arg_hash->{resultsdir} = "resultsAnneal";
$arg_hash->{numannealedruns} = 100;
#$arg_hash->{samples}  = 6000;
#$arg_hash->{burnin}   = 1000;
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
    unless (-e "$args->{resultsdir}"){
	mkdir("$args->{resultsdir}");
    }
    system($command);

# Comment out the next three lines to run admixmap without R script
    print "Starting R script to process output\n";
    system("R --quiet --no-save --no-restore <AdmixmapOutput.R >results/Rlog.txt RESULTSDIR=$args->{resultsdir}");
    print "R script completed\n\n";
}



