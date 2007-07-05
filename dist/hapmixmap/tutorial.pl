#!/usr/bin/perl
use strict; 
use File::Path;
use Getopt::Long;

my $function_file = "doanalysis.pl";

my $train = '';
my $test = '';
my $resume = '';
my $usage = 0;
my $executable = "hapmixmap";
my $rscript = "AdmixmapOutput.R";

GetOptions("train" =>\$train,
	   "test" => \$test,
	   "resume" => \$resume,
	   "help!" => \$usage,
           "function-file=s" => \$function_file,
           "exec=s"          => \$executable,
           "rscript=s"       => \$rscript
	  );

if (not ($train or $test) or $usage) {
    print "Usage: $0 < --train | --test > [--resume] [--exec=...] [--rscript=...] [--function-file=...]\n";
}

##the following lines are a botch to make the script work straight out of the repository
if(!(-f $function_file) && (-f "../$function_file")){
##try one level up
  $function_file = "../$function_file";
}
if(!(-f $rscript) && (-f "../$rscript")){
##try one level up
  $rscript = "../$rscript";
}

require $function_file or die("cannot find doanalysis.pl");

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
		arrivalrateprior => "12, 0.1, 2000,1000",
		residualadprior => "0.25, 1",
		arrivalratesamplerparams => "0.1, 0.00001, 10, 0.9, 20",
		
		resultsdir => "results",
		logfile => "logfile.txt",
		paramfile => "paramfile.txt",
		residualadfile => "allelefreqdispersion.txt",
		ergodicaveragefile => "ergodicaverage.txt"
};

if($train){
  # Initial run with HapMap data
  $arg_hash->{resultsdir} = "ResultsTraining";
  $arg_hash->{finalvaluedir} = "data";
  if($resume){#use initial values from previous run
    $arg_hash->{initialvaluedir} = "data";
  }
  doAnalysis($executable, $rscript, $arg_hash);
}

if($test){
  $arg_hash->{finalvaluedir} = "initialValues";
  if($resume){#use initial values from previous run
    $arg_hash->{initialvaluedir} = "initialValues";
  }else{#use initial values from training run
    $arg_hash->{initialvaluedir} = "data";
  }
  # Run with HapMap + Case/Control data
  $arg_hash->{resultsdir} = "Results";
  doAnalysis($executable, $rscript, $arg_hash);
}
