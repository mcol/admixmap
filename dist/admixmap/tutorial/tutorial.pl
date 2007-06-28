#!/usr/bin/perl
use strict; 
use File::Path;

my $function_file = "../doanalysis.pl";

require $function_file or die("cannot find doanalysis.pl");

# Change this to the location of the admixmap executable
my $executable = '../admixmap';

# Change this to the location of the R script
my $rscript = "../AdmixmapOutput.R";

my $arg_hash = {
#data files
    genotypesfile                   => 'data/genotypes.txt',
    locusfile                       => 'data/loci.txt',
    covariatesfile                  => 'data/covariates2std.txt', # age, sex 
    outcomevarfile                  => 'data/outcomevars.txt',
#main options
    samples  => 1200,
    burnin   => 200,
    every    => 5,
    numannealedruns => 0, #200, # 100, 
    displaylevel => 2,
#output file options
    logfile                     => 'log.txt',
    regparamfile                => 'regparam.txt',
    ergodicaveragefile          => 'cumulativeAverages.txt',
# optional tests
    residualldtest               => 1,
    allelicassociationtest       => 1,
    haplotypeassociationtest     => 1,
    stratificationtest           => 1,
    hwtest                       => 1
};

# model with reference prior on allele freqs in 1 population, skin reflectance as continuous outcome var
$arg_hash->{populations}           = 1;
$arg_hash->{resultsdir}            = 'SinglePopResults';
$arg_hash->{outcomevarcols}       = 2; # skin reflectance
#&doAnalysis($executable, $rscript, $arg_hash);

# model with reference prior on allele freqs in 2 populations
$arg_hash->{populations}           = 2;
$arg_hash->{samples}   = 6000;
$arg_hash->{burnin}    = 1000;
$arg_hash->{paramfile}                 = 'popadmixparams.txt',
$arg_hash->{resultsdir}            = 'TwoPopsResults';  
#&doAnalysis($executable, $rscript, $arg_hash);

# model with reference prior on allele freqs in 3 populations
$arg_hash->{populations}           = 3;
$arg_hash->{resultsdir}            = 'ThreePopsResults';  
#&doAnalysis($executable, $rscript, $arg_hash);

# model with prior allele freqs 
delete $arg_hash->{populations};
$arg_hash->{resultsdir}                    = 'PriorFreqResultsSkin';  
$arg_hash->{samples}    = 1200;
$arg_hash->{burnin}    = 200;
$arg_hash->{priorallelefreqfile}           = 'data/priorallelefreqs.txt';
$arg_hash->{dispersiontest}            = 1;
$arg_hash->{indadmixturefile}              = 'indivadmixture.txt';
$arg_hash->{ancestryassociationtest}  = 1;
&doAnalysis($executable, $rscript, $arg_hash);

# model with prior allele freqs and diabetes as binary outcome var 
delete $arg_hash->{populations};
$arg_hash->{resultsdir}                = 'PriorFreqResultsDiabetes';  
$arg_hash->{outcomevarcols}           = 1; # diabetes as outcome
$arg_hash->{affectedsonlytest}    = 1
#$arg_hash->{thermo}    = 1;
#$arg_hash->{numannealedruns}    = 100;
&doAnalysis($executable, $rscript, $arg_hash);

# model with fixed allele freqs and diabetes as binary outcome var 
delete $arg_hash->{populations};
$arg_hash->{resultsdir}                = 'FixedAlleleFreqResultsDiabetes';  
$arg_hash->{affectedsonlytest}   = 1;
$arg_hash->{fixedallelefreqs}    = 1;
#$arg_hash->{numannealedruns}    = 100;
doAnalysis($executable, $rscript, $arg_hash);

# model with historic allele freqs and both outcome vars
delete $arg_hash->{targetindicator};  
delete $arg_hash->{priorallelefreqfile};
delete $arg_hash->{dispersiontestfile};
delete $arg_hash->{affectedsonlytest};
$arg_hash->{fixedallelefreqs}    = 0;
$arg_hash->{numannealedruns}    = 0;
$arg_hash->{resultsdir}                = 'HistoricFreqResults';  
$arg_hash->{randommatingmodel}         = 1;
$arg_hash->{outcomevarcols}           = "1,2"; # both outcome vars
$arg_hash->{historicallelefreqfile}    = "data/priorallelefreqs.txt";
$arg_hash->{etapriorfile}              = "data/etapriors.txt";
$arg_hash->{dispparamfile}             = "dispersionparams.txt";
$arg_hash->{fstoutputfile}             = "lociFst.txt";
$arg_hash->{allelefreqoutputfile}      = "allelefreqs.txt";
&doAnalysis($executable, $rscript, $arg_hash);

