#!/usr/bin/perl -w
use strict;
use File::Path;

my $DEBUG = 0; # zero gives less output

# Change this to the location of the admixmap executable
my $executable = 'c:\cvs\genepi\test\admixmap3.2b';

# $arg_hash is a hash of parameters passed to
# the executable as arguments.
#
# keys (left-hand side) are parameter names
# values (right-hand side) are parameter values
my $arg_hash = 
{
#data files
    genotypesfile                   => 'data/genotypes.txt',#no missing data
    #genotypesfile                   => 'regdata/simgenotypes2.txt', #with missing data

    locusfile                          => 'data/loci.txt',
    #priorallelefreqfile             => 'data/priorallelefreqs.txt',
 priorallelefreqfile =>'data/trueallelefreqs1pop.txt',
    #fixedallelefreqs => 1,
    globalrho => 1,
#populations =>2,
    #randommatingmodel =>1,

    #covariatesfile                  => 'regdata/covariates3.txt',
    #outcomevarfile               => 'regdata/binoutcome.txt',#binary outcome
    #outcomevarfile               => 'data/binoutcome.txt', #continuous outcome

#main options
    displaylevel   => 3, #verbose output
    #targetindicator => 0, # column 1
    samples  => 250,
    burnin   => 50,
    every    => 1,
thermo => 0,
numannealedruns => 0,


#output files
    resultsdir               => "resultsoldexec1pop",
    logfile                     => 'logfile.txt',
    paramfile               => 'paramfile.txt',
    regparamfile          => 'regparamfile.txt',
    indadmixturefile     => 'indadmixture.txt',
    ergodicaveragefile => 'ergodicaverage.txt',
    allelefreqoutputfile  => 'allelefreqoutputfile.txt',

#optional tests
 #residualallelicassocscorefile => 'residualallelicscoretests.txt',
    #hwtest => 'HWTest.txt',
#dispersiontestfile =>'dispersiontest.txt',
    #allelicassociationscorefile       => 'allelicassociationscorefile.txt',
    #ancestryassociationscorefile  => 'ancestryassociationscorefile.txt',
    #affectedsonlyscorefile             => 'affectedsonlyscorefile.txt',
    #haplotypeassociationscorefile => 'hapassocscore.txt',
    #stratificationtestfile                   => 'strat_test.txt'
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

    $ENV{'RESULTSDIR'} = $args->{resultsdir};
    print "Results will be written to subdirectory $ENV{'RESULTSDIR'}\n";
    system($command);
    print "Starting R script to process output\n";
    system("R CMD BATCH --quiet --no-save --no-restore c:/cvs/genepi/test/AdmixmapOutput.R $args->{resultsdir}//Rlog.txt");
    #system("R CMD BATCH --quiet --no-save --no-restore /home/david/genepi/test/AdmixmapOutput.R $args->{resultsdir}//Rlog.txt");
    print "R script completed\n\n";
}



