#!/usr/bin/perl
use strict; 
use File::Path;
use Getopt::Long;

my $function_file = "doanalysis.pl";

require $function_file or die("cannot find doanalysis.pl");

my $train = '';
my $test = '';
my $resume = '';

GetOptions("train" =>\$train,
	   "test" => \$test,
	   "resume" => \$resume,
	  );


my $executable = "./hapmixmap";
my $rscript = "./AdmixmapOutput.R";
################### DO NOT EDIT ABOVE THIS LINE ########################


my $arg_hash = {
genotypesfile => "data/training_genotypes.txt",
locusfile => "data/training_loci.txt",
checkdata => 0,

samples => 50,
burnin => 10,
every => 1,
numannealedruns => 0,
thermo => 0,
states => 8,
hapmixmodel =>  1,
displaylevel => 3,

freqdispersionhiermodel => 1,
allelefreqprior => "0.2, 1, 1",
hapmixlambdaprior => "150, 1, 40, 10",
lambdasamplerparams => "0.1, 0.00001, 10, 0.9, 20",

resultsdir => "results",
logfile => "logfile.txt",
paramfile => "paramfile.txt",
dispparamfile => "allelefreqpriors.txt",

#allelefreqprioroutputfile => "initialetas.txt",
#allelefreqoutputfile => "initialallelefreqs.txt",
#hapmixlambdaoutputfile => "initiallambdas.txt",

ergodicaveragefile => "ergodicaverage.txt"
};

if($train){
  # Initial run with HapMap data
  $arg_hash->{resultsdir} = "ResultsTraining";
  if($resume){#use initial values from previous run
    $arg_hash->{initialallelefreqfile} = "ResultsTraining/state-allelefreqs.txt";
    $arg_hash->{initialfreqpriorfile} = "ResultsTraining/state-freqpriors.txt";
    $arg_hash->{initiallambdafile} => "ResultsTraining/state-lambdas.txt";
  }
  doAnalysis($executable, $rscript, $arg_hash);
}

if($test){
  if($resume){#use initial values from previous run
    $arg_hash->{initialallelefreqfile} = "Results/state-allelefreqs.txt";
    $arg_hash->{initialfreqpriorfile} = "Results/state-freqpriors.txt";
    $arg_hash->{initiallambdafile} => "Results/state-lambdas.txt";
  }else{#use initial values from traiing run
    $arg_hash->{initialallelefreqfile} = "ResultsTraining/state-allelefreqs.txt";
    $arg_hash->{initialfreqpriorfile} = "ResultsTraining/state-freqpriors.txt";
    $arg_hash->{initiallambdafile} => "ResultsTraining/state-lambdas.txt";
  }
  # Run with HapMap + Case/Control data
  $arg_hash->{resultsdir} = "Results";
  doAnalysis($executable, $rscript, $arg_hash);
}
