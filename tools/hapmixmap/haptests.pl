#!/usr/bin/perl
use strict; 
use File::Path;
use Getopt::Long;

my $parallel = '';
my $simulate = '';
my $samples  = 25;
my $burnin = 5;
my $every = 1;

GetOptions("parallel" =>\$parallel, 
	   "simulate" => \$simulate, 
	   "samples=i" => \$samples,
	   "burnin=i" => \$burnin,
	   "every=i" => \$every);

if($simulate){ 
    print "Running R script to simulate data\n";
    system("R CMD BATCH --vanilla simHapMix.R");
    print "simulation complete\n";
}

my $function_file = "../doanalysis.pl";

require $function_file or die("cannot find doanalysis.pl");

my $serial_executable = "$ENV{'HOME'}/bin/hapmixmap";
my $parallel_executable = "mpiexec $ENV{'HOME'}/bin/hapmixmap-para";
my $rscript = "../admixmap/AdmixmapOutput.R";
################### DO NOT EDIT ABOVE THIS LINE ########################


my $executable = $serial_executable;
if($parallel){
    $executable = $parallel_executable;
}

my $arg_hash = {
#data files
    locusfile                       => 'data/loci.txt',

    #genotypesfile                   => 'data/genotypes.txt', #diploid data
     genotypesfile                   => 'data/genotypes_haploid.txt',#haploid data

    #priorallelefreqfile             => 'data/allelefreqs.txt',
    #fixedallelefreqs                => 1,

    states=>4,

#main options
    resultsdir => 'results',
    displaylevel   => 3, 

    samples  => $samples,
    burnin   => $burnin,
    every    => $every,

    numannealedruns => 0,
    thermo => 0,
    hapmixmodel => 1,

    checkdata=> 0,

hapmixlambdaprior=>"400, 1, 10, 1",

allelefreqprior => "2, 10, 1",
freqdispersionhiermodel => 1,

#initialhapmixlambdafile => "data/initialambdas.txt",
#allelefreqfile => "data/initialallelefreqs.txt",

lambdasamplerparams => "0.5, 0.00001, 10, 0.9, 20",

#output files
    logfile                     => 'logfile.txt',
    paramfile               => 'paramfile.txt',
    dispparamfile => "allelefreqpriorsamples.txt",
    #regparamfile          => 'regparamfile.txt',
    ergodicaveragefile => 'ergodicaverage.txt',

#final values
    finalallelefreqfile  => "initialallelefreqs.txt",
    finalfreqpriorfile =>"initialallelefreqpriors.txt",
    finallambdafile =>"initiallambdas.txt",

#posterior means
    hapmixlambdaoutputfile => "lambdaPosteriorMeans.txt",
    allelefreqprioroutputfile => "freqDispersionPosteriorMeans.txt"

#optional tests
#residualallelicassocscorefile => 'residualLDscores.txt',
    #allelicassociationscorefile       => 'allelicassociationscorefile.txt',
};

# Initial run 
#haploid data
$arg_hash->{resultsdir}            = 'Results_Haploid';
doAnalysis($executable, $rscript, $arg_hash);

#diploid data
#$arg_hash->{resultsdir}            = 'ResultsDiploid';
#$arg_hash->{genotypesfile} = "data/genotypes.txt";
#doAnalysis($executable,$arg_hash);
#system("cp Results/initialallelefreqs.txt data");

# rerun with final values of previous run as intial values of this
#system("cp $arg_hash->{resultsdir}/initialallelefreqs.txt data");
$arg_hash->{resultsdir}            = 'ResultsRerun';
$arg_hash->{initialallelefreqfile}="data/initialallelefreqs.txt";
$arg_hash->{initialhapmixlambdafile}="data/initiallambdas.txt";
$arg_hash->{initialfreqpriorfile} = "data/initialfreqpriors.txt";
#$arg_hash->{fixedallelefreqs} = 0;
#delete $arg_hash->{priorallelefreqfile};
#doAnalysis($executable,$arg_hash);
