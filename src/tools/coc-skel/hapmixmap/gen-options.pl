#!/usr/bin/perl
# script to generate options files for multiple runs of hapmixmap with different prior or different numbers of states

sub writeOptionsFile{
    my $hash = $_[0];
    my $filename = $_[1];
    open(OPTIONFILE, ">$filename") or die ("Could not open options file");
    foreach my $key (keys %$hash){
	print OPTIONFILE $key . '=' . $hash->{$key} . "\n";
    }
    close OPTIONFILE;
};


###################################################################
my $dataprefix="data/chr22/hapmixmap";
my @Panels=("CEU", "YRI", "JPTCHB");
my $resultsprefix="results2";

my $arg_hash = {
    deleteoldresults => 0,

#model
    states                => 6,
    checkdata             => 0,
    fixedmixtureprops     => 1,
    fixedmixturepropsprecision => 1,

    #main options
    displaylevel    => 3,
    samples         => 12, #$samples,
    burnin          => 2,  #$burnin,
    every           => 1,  #$every,
    numannealedruns => 0,
    thermo          => 0,
    hapmixmodel     => 1,

    #prior spec
    freqprecisionhiermodel => 0,
    arrivalratesamplerparams => "0.1, 0.00001, 10, 0.9, 20",
    
    #output files
    logfile =>'logfile.txt',
    paramfile         =>'paramfile.txt',#mean and var of sampled arrival rates
    freqprecisionfile =>'freqprecision.txt', #mean and var of sampled allele freq precision
    arrivalrateposteriormeanfile => "ArrivalRatePosteriorMeans.txt",
};


for my $pop(@Panels) { # loop over 3 populations
    print "$pop\n";
    
    # data files
    $arg_hash->{locusfile}="$dataprefix/$pop/phased_5000_loci.txt";
    $arg_hash->{genotypesfile}="$dataprefix/$pop/phased_5000_genotypes.txt";
    
    ## priors
    $arg_hash->{arrivalrateprior}="1.2, 0.05, 1, 0.5";
    $arg_hash->{residualadprior}="2, 8";
    
    $run_name = 
	"$pop"."\_"."arp-$arg_hash->{arrivalrateprior}_afpp-$arg_hash->{residualadprior}";
    $run_name =~ s/ //g;
    $run_name =~ s/\,/\-/g;
    print "$run_name\n";
    $arg_hash->{resultsdir}="$resultsprefix/$run_name";
    $arg_hash->{finalvaluedir}="$resultsprefix/$run_namevalues";
    
    # training run
    $arg_hash->{predictgenotypes} = 0;
    my $optionsfilename = "training_$run_name.txt";
    writeOptionsFile($arg_hash, $optionsfilename);
    system("hapmixmap $optionsfilename");
    
    # testing run
    $arg_hash->{testgenotypesfile}="$dataprefix/$pop/obs_masked_genotypes.txt";
    $optionsfilename = "testing_$run_name.txt";
    $arg_hash->{predictgenotypes} = 1;
    writeOptionsFile($arg_hash, $optionsfilename);
    system("hapmixmap $optionsfilename");
}
