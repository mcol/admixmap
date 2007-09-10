#!/usr/bin/perl
# script to generate options files for multiple runs of hapmixmap with different prior or different numbers of states
use strict; 

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
    # training run
    $arg_hash->{predictgenotypes} = 0;
    if(exists $arg_hash->{testgenotypesfile}) {
	delete $arg_hash->{testgenotypesfile};
    }
    my $optionsfilename = "training_$run_name.txt";
    writeOptionsFile($arg_hash, $optionsfilename);
    system("hapmixmap $optionsfilename");
    
    # testing run
    $arg_hash->{testgenotypesfile}="$testgenotypesfile";
    $optionsfilename = "testing_$run_name.txt";
    $arg_hash->{predictgenotypes} = 1;
    writeOptionsFile($arg_hash, $optionsfilename);
    system("hapmixmap $optionsfilename");
};
  
###################################################################
my $dataprefix="data/chr22/hapmixmap";
my @Panels=("YRI", "CEU", "JPTCHB");
my $resultsprefix="results2";

my $arg_hash = {
    deleteoldresults => 0,

#model
    states                => 6,
    checkdata             => 0,
    fixedmixtureprops     => 1,
    fixedmixturepropsprecision => 1,

    #main options
    displaylevel    => 2,
    samples         => 250, #$samples,
    burnin          => 50,  #$burnin,
    every           => 5,   #$every,
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
    
    # data files
    $arg_hash->{locusfile}="$dataprefix/$pop/train_loci.txt";
    $arg_hash->{genotypesfile}="$dataprefix/$pop/train_genotypes.txt";
    my $testgenotypesfile="$dataprefix/$pop/test_genotypes.txt";
    
    $arg_hash->{arrivalrateprior} = "1.2, 0.05, 1, 0.5";
    $arg_hash->{residualadprior} = "0.25, 1";

    # loop over numbers of states
    my @numstates = (2, 4, 6, 8, 10);
    for my $num(@numstates) {
        $arg_hash->{states} = $num;
	my $run_name = "$pop"."\_".
	    "states-$arg_hash->{states}"."\_".
	    "arp-$arg_hash->{arrivalrateprior}"."\_".
	    "afpp-$arg_hash->{residualadprior}";
	$run_name =~ s/ //g;
	$run_name =~ s/\,/\-/g;
	print "$run_name\n";
	
	# output folders
	$arg_hash->{resultsdir}="$resultsprefix/$run_name";
	$arg_hash->{finalvaluedir}="$resultsprefix/$run_name"."_fv";
	
	# write config files
	trainandtest($arg_hash, $run_name, $testgenotypesfile);
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
		= "$pop"."\_".
		"states-$arg_hash->{states}"."\_".
		"arp-$arg_hash->{arrivalrateprior}"."\_".
		"afpp-$arg_hash->{residualadprior}";
	    $run_name =~ s/ //g;
	    $run_name =~ s/\,/\-/g;
	    print "$run_name\n";
	    
	    # output folders
	    $arg_hash->{resultsdir}="$resultsprefix/$run_name";
	    $arg_hash->{finalvaluedir}="$resultsprefix/$run_name"."_fv";
	    
	    # write config files
	    trainandtest($arg_hash, $run_name, $testgenotypesfile);
	}
    }
};
