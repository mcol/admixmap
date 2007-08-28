# script to generate options files for multiple runs of hapmixmap with different prior or different numbers of states

my $arg_hash = {
    deleteoldresults => 0,

#data
    genotypesfile                   => "$datadir/phased_genotypes.txt",
    locusfile                       => "$datadir/phased_loci.txt",

#model
    states     => 6,
    checkdata       => 0,
    # mixturepropsprior => "50, 1",
    fixedmixtureprops => 1,
    fixedmixturepropsprecision =>1,

#main options
    resultsdir      => 'results',
    displaylevel    => 3,
    samples         => $samples,
    burnin          => $burnin,
    every           => $every,
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

## TO BE CONTINUED ...
