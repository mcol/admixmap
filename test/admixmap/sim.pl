#!/usr/bin/perl -w
use strict;

my $function_file = "../../dist/doanalysis.pl";#"$ENV{'HOME'}/genepi/trunk/dist/doanalysis.pl";
require $function_file or die("cannot find doanalysis.pl");

my $executable = 'admixmap';
my $rscript = "../../dist/AdmixmapOutput.R";

my $arg_hash = {
#data files
    genotypesfile                   => 'simdata/genotypes.txt',
    locusfile                       => 'simdata/loci.txt',
    priorallelefreqfile             => 'simdata/priorallelefreqs.txt',
    # covariatesfile                  => 'data/covariates3.txt',
    outcomevarfile                  => 'simdata/outcome.txt',
    #fixedallelefreqs                => 1,
#main options
    displaylevel   => 3, #verbose output
    samples  => 250,
    burnin   => 50,
    every    => 5,
    populations => 2,
    numannealedruns => 0,
    #indadmixhiermodel => 0,
    #admixtureprior => "3,1",

#output files
    resultsdir               => "simResults",
    logfile                  => 'logfile.txt',
    paramfile                => 'paramfile.txt',
    regparamfile             => 'regparam.txt',
    indadmixturefile         => 'indadmixture.txt',
    #ergodicaveragefile       => 'ergodicaverage.txt',
    allelefreqoutputfile  => 'allelefreqoutputfile.txt',

#optional tests
    allelicassociationtest   => 1,
    ancestryassociationtest  => 1,
    #affectedsonlytest        => 1,
    #haplotypeassociationtest => 1,
    dispersiontest           => 1,
    hwtest                   => 1,
    #stratificationtest       => 1,
    residualldtest            => 1
};

doAnalysis($executable, $rscript, $arg_hash);

