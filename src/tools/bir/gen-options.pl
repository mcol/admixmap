#!/usr/bin/perl
# script to generate options files for multiple runs of hapmixmap with different prior or different numbers of states
use strict; 
use File::Path;

sub writeOptionsFile{
    my $hash = $_[0];
    my $filename = $_[1];
    open(OPTIONFILE, ">$filename") or die ("Could not open options file");
    foreach my $key (keys %$hash){
	print OPTIONFILE $key . '=' . $hash->{$key} . "\n";
    }
    close OPTIONFILE;
};

sub trainandtest {
    my $arg_hash = $_[0];
    my $run_name = $_[1];
    my $testgenotypesfile = $_[2];
    my $comp_states = $_[3]; #1 if comparing states
    # training run
    print "training run\n";
    if(exists $arg_hash->{initialvaluedir}) {
        delete $arg_hash->{initialvaluedir};
    }
    $arg_hash->{predictgenotypes} = 0;
    if(exists $arg_hash->{testgenotypesfile}) {
	delete $arg_hash->{testgenotypesfile};
    }
    my $optionsfilename = "configfiles/training_$run_name.txt";
    writeOptionsFile($arg_hash, $optionsfilename);
    if($comp_states ==1){
	print TRAIN_COMP_STATES_LIST "hapmixmap $optionsfilename\n";
    }else{
	print TRAIN_COMP_PRIORS_LIST "hapmixmap $optionsfilename\n";
    }
    system("hapmixmap $optionsfilename");
    
    # testing run
    print "testing run\n";
    $arg_hash->{testgenotypesfile}="$testgenotypesfile";
    $arg_hash->{initialvaluedir}=$arg_hash->{finalvaluedir};
    $arg_hash->{predictgenotypes} = 1;
    $optionsfilename = "configfiles/testing_$run_name.txt";
    writeOptionsFile($arg_hash, $optionsfilename);
    if($comp_states ==1){
	print TEST_COMP_STATES_LIST "hapmixmap $optionsfilename\n";
    }else{
	print TEST_COMP_PRIORS_LIST "hapmixmap $optionsfilename\n";
    }
    system("hapmixmap $optionsfilename");
};

###################################################################

## change these prefixes to run entire chromosome
my $whichchr = "chr22_5kloci";
#my $whichchr = "chr22"; 

my $dataprefix="data/$whichchr/hapmixmap"; 
my $resultsprefix="results/$whichchr/hapmixmap"; 

my @Panels=("YRI", "CEU", "JPTCHB");
my @seeds=(2190, 3367, 5211, 7318);

if(!(-e "configfiles")) { 
    mkpath "configfiles" or die "cannot make directory configfiles";
}

open(TRAIN_COMP_STATES_LIST, ">compare_states_train_tasks.txt") or die ("could not open task list");
open(TEST_COMP_STATES_LIST, ">compare_states_test_tasks.txt") or die ("could not open task list");
open(TRAIN_COMP_PRIORS_LIST, ">compare_priors_train_tasks.txt") or die ("could not open task list");
open(TEST_COMP_PRIORS_LIST, ">compare_priors_test_tasks.txt") or die ("could not open task list");


my $arg_hash = {
    deleteoldresults => 1,

#model
    states                => 6,
    checkdata             => 0,
    fixedmixtureprops     => 1,
    fixedmixturepropsprecision => 1,

    #main options
    displaylevel    => 0,
    samples         => 12, #$samples,
    burnin          => 2,  #$burnin,
    every           => 1,   #$every,
    numannealedruns => 0,
    thermo          => 0,
    hapmixmodel     => 1,

    #prior spec
    freqprecisionhiermodel => 0,
    arrivalratesamplerparams => "0.1, 0.00001, 10, 0.9, 20",
    
    #output files
    logfile =>'logfile.txt',
    paramfile         =>'paramfile.txt', #mean and var of sampled arrival rates
    freqprecisionfile =>'freqprecision.txt', #mean and var of sampled allele freq precision
    arrivalrateposteriormeanfile => "ArrivalRatePosteriorMeans.txt",
};

print "\n";
for my $pop(@Panels) { # loop over 3 populations
  my $s = 1;
  for my $seed(@seeds[0]){
  #for my $seed(@seeds){
    # data files
    $arg_hash->{locusfile}="$dataprefix/$pop/train_loci.txt";
    $arg_hash->{genotypesfile}="$dataprefix/$pop/train_genotypes.txt";
    $arg_hash->{seed}=$seed;
    my $testgenotypesfile="$dataprefix/$pop/test_genotypes.txt";
    
    $arg_hash->{arrivalrateprior} = "1.2, 0.05, 1, 0.5";
    $arg_hash->{residualadprior} = "0.25, 1";
    
    # loop over numbers of states
    my @numstates = (2, 4, 6, 8, 10);
    for my $num(@numstates) {
      $arg_hash->{states} = $num;
      my $run_name = 
	"states-$arg_hash->{states}"."\_".
	  "arp-$arg_hash->{arrivalrateprior}"."\_".
	    "afpp-$arg_hash->{residualadprior}"."\_".
	      "seed-$s";;
      $run_name =~ s/ //g;
      $run_name =~ s/\,/\-/g;
      print "$run_name\n";
      
      # output folders
      $arg_hash->{resultsdir}="$resultsprefix/$pop/$run_name";
      $arg_hash->{finalvaluedir}="$resultsprefix/$pop/$run_name"."_fv";
      
      # write config files
      trainandtest($arg_hash, $pop."_".$run_name, $testgenotypesfile, 1);
    }
    
    $arg_hash->{states} = 6; # default
    # loop over priors
    my @apriors = ("1.2, 0.05, 1, 0.5", "1.2, 0.05, 2", 
		   "12, 0.5, 1, 0.5",   "12, 0.5, 2");
    for my $aprior(@apriors) { # loop over arrival rate priors
      
      $arg_hash->{arrivalrateprior}= $aprior;
      
      my @rpriors = ("0.25, 1", "2, 8");
      for my $rprior(@rpriors) { # loop over allelic diversity priors
	$arg_hash->{residualadprior}=$rprior;
	
	my $run_name 
	  = "states-$arg_hash->{states}"."\_".
	    "arp-$arg_hash->{arrivalrateprior}"."\_".
	      "afpp-$arg_hash->{residualadprior}"."\_".
		"seed-$s";
	$run_name =~ s/ //g;
	$run_name =~ s/\,/\-/g;
	print "$run_name\n";
	
	# output folders
	$arg_hash->{resultsdir}="$resultsprefix/$pop/$run_name";
        if( !(-e "$arg_hash->{resultsdir}") ) {
	    mkpath "$arg_hash->{resultsdir}" 
		or die "cannot make directory $arg_hash->{resultsdir}\n"; 
	};
	$arg_hash->{finalvaluedir}="$resultsprefix/$pop/$run_name"."_fv";
	
	# write config files
	trainandtest($arg_hash, $pop."_".$run_name, $testgenotypesfile, 0);
      }
    }
    $s = $s+1;
  }
};

close(TRAIN_COMP_STATES_LIST);
close(TEST_COMP_STATES_LIST);
close(TRAIN_COMP_PRIORS_LIST);
close(TEST_COMP_PRIORS_LIST);
