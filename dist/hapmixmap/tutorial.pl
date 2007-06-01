#!/usr/bin/perl
use strict; 
use File::Path;
use Getopt::Long;

my $function_file = "./doanalysis.pl";

require $function_file or die("cannot find doanalysis.pl");

my $train = '';
my $test = '';
my $resume = '';
my $usage = 0;

GetOptions("train" =>\$train,
	   "test" => \$test,
	   "resume" => \$resume,
	   "help!" => \$usage,
	  );

if (not ($train or $resume or $test) or $usage) {
    print "Usage: $0 [ --train | --resume | --test ]\n";
}


my $executable = "./hapmixmap";
my $rscript = "./AdmixmapOutput.R";


my $arg_hash = {
genotypesfile => "data/dn_genotypes.txt",
locusfile => "data/dn_loci.txt",
checkdata => 0,

samples => 50,
burnin => 10,
every => 1,
numannealedruns => 0,
thermo => 0,
states => 8,
hapmixmodel =>  1,
displaylevel => 3,

residualadhiermodel => 0,
arrivalrateprior => 12. 0.1, 2000,1000,
residualadprior => 0.25, 1,
arrivalratesamplerparams => "0.1, 0.00001, 10, 0.9, 20",

resultsdir => "results",
logfile => "logfile.txt",
paramfile => "paramfile.txt",
residualadfile => "allelefreqdispersion.txt",

#allelefreqprioroutputfile => "initialetas.txt",
#allelefreqoutputfile => "initialallelefreqs.txt",
#arrivalrateoutputfile => "initiallambdas.txt",

ergodicaveragefile => "ergodicaverage.txt"
};

if($train){
  # Initial run with HapMap data
  $arg_hash->{resultsdir} = "ResultsTraining";
  if($resume){#use initial values from previous run
    $arg_hash->{initialallelefreqfile} = "ResultsTraining/state-allelefreqs.txt";
    $arg_hash->{initialfreqpriorfile} = "ResultsTraining/state-freqpriors.txt";
    $arg_hash->{initialarrivalratefile} => "ResultsTraining/state-arrivalrates.txt";
  }
  doAnalysis($executable, $rscript, $arg_hash);
}

if($test){
  if($resume){#use initial values from previous run
    $arg_hash->{initialallelefreqfile} = "Results/state-allelefreqs.txt";
    $arg_hash->{initialfreqpriorfile} = "Results/state-freqpriors.txt";
    $arg_hash->{initialarrivalratefile} => "Results/state-arrivalrates.txt";
  }else{#use initial values from traiing run
    $arg_hash->{initialallelefreqfile} = "ResultsTraining/state-allelefreqs.txt";
    $arg_hash->{initialfreqpriorfile} = "ResultsTraining/state-freqpriors.txt";
    $arg_hash->{initialarrivalratefile} => "ResultsTraining/state-arrivalrates.txt";
  }
  # Run with HapMap + Case/Control data
  $arg_hash->{resultsdir} = "Results";
  doAnalysis($executable, $rscript, $arg_hash);
}
