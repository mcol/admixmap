#!/usr/bin/perl -w
use strict;

my $DEBUG = 0; # zero gives less output

# Change this to the location of the admixmap executable
my $executable = './admixmap';
# Change this to desired name of results directory
my $resultsdir = "CattleResults";

# $arg_hash is a hash of parameters passed to
# the executable as arguments.
#
# keys (left-hand side) are parameter names
# values (right-hand side) are parameter values
my $arg_hash = 
{
#data files
    genotypesfile                   => 'Cattledata/ngu_complete2.txt',
    locusfile                          => 'Cattledata/loci.txt',
    populations => 2,

#main options
 globalrho => 1,
#randommatingmodel => 0,

    samples  => 2500,
    burnin   => 500,
    every    => 50,

#output files
    resultsdir               => "$resultsdir",
    logfile                     => 'log.txt',
    paramfile               => 'paramfile.txt',
    indadmixturefile     => 'indadmixture.txt',
    ergodicaveragefile => 'ergodicaverages.txt',
    allelefreqoutputfile  => 'allelefreqs.txt',

#optional tests
#admixturescorefile => 'admixscorefile.txt',
    #allelicassociationtest   => 1,
    #ancestryassociationtest  => 1,
    #affectedsonlytest        => 1,
    #haplotypeassociationtest => 1,
    #stratificationtest       => 1
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



