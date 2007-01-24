#!/usr/bin/perl
use strict; 
use File::Path;
use Getopt::Long;

my $parallel = '';
my $simulate = '';

GetOptions("parallel" =>\$parallel, "simulate" => \$simulate);

sub getArguments
{
    my $hash = $_[0];
#    my $arg = '';
#   foreach my $key (keys %$hash){
#	$arg .= ' --'. $key .'='. $hash->{$key};
#    }
#    return $arg;
    my $filename = 'perlargs.txt';
    open(OPTIONFILE, ">$filename") or die ("Could not open args file");
    foreach my $key (keys %$hash){
      print OPTIONFILE $key . '=' . $hash->{$key} . "\n";
    }
    close OPTIONFILE;
    return " ".$filename;
}

sub doAnalysis
{
    if($simulate){ 
	print "Running R script to simulate data\n";
	system("R CMD BATCH --vanilla simHapMix.R");
	print "simulation complete\n";
    }
    my ($prog,$args) = @_;
    if($parallel){
	$args->{resultsdir} = "$args->resultsdir". "Parallel";
    }
    my $command = "";
    $parallel ? $command = "mpiexec " : $command =  "";
    $command = $command . $prog.getArguments($args);
    
    $ENV{'RESULTSDIR'} = $args->{resultsdir};
    print "\nResults will be written to subdirectory $ENV{'RESULTSDIR'}\n";
    print $command;
    if(system($command) ==0){
	print "Starting R script to process output\n";
	system("R CMD BATCH --quiet --no-save --no-restore ../test/AdmixmapOutput.R $args->{resultsdir}/Rlog.txt");
	print "R script completed\n\n";
    }
}

################### DO NOT EDIT ABOVE THIS LINE ########################
my $serial_executable = '../test/hapmixmap';
my $parallel_executable = '../test/hapmixmap-para';
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

    samples  => 25,
    burnin   => 5,
    every    => 1,

    numannealedruns => 0,
    thermo => 0,
    hapmixmodel => 1,
#   indadmixhiermodel => 0,
    randommatingmodel => 0,
    checkdata=> 0,

hapmixlambdaprior=>"400, 1, 10, 1",

allelefreqprior => "2, 10, 1",
#freqdispersionhiermodel => 0,

#initialhapmixlambdafile => "data/initialambdas.txt",
#allelefreqfile => "data/initialallelefreqs.txt",

rhosamplerparams => "0.5, 0.00001, 10, 0.9, 20",

#output files
    logfile                     => 'logfile.txt',
    paramfile               => 'paramfile.txt',
    dispparamfile => "alelefreqpriorsamples.txt",
    #regparamfile          => 'regparamfile.txt',
    ergodicaveragefile => 'ergodicaverage.txt',

    allelefreqoutputfile  => "initialallelefreqs.txt",
    allelefreqprioroutputfile =>"allelefreqpriors.txt",
    hapmixlambdaoutputfile => "data/initiallambdas.txt",

#optional tests
#residualallelicassocscorefile => 'residualLDscores.txt',
    #allelicassociationscorefile       => 'allelicassociationscorefile.txt',
};

# Initial run 
#haploid data
$arg_hash->{resultsdir}            = 'ResultsHaploid';  
doAnalysis($executable,$arg_hash);

#diploid data
$arg_hash->{resultsdir}            = 'ResultsDiploid';  
$arg_hash->{genotypesfile} = "data/genotypes.txt";
#doAnalysis($executable,$arg_hash);
#system("cp Results/initialallelefreqs.txt data");

# rerun with final values of previous run as intial values of this
#system("cp $arg_hash->{resultsdir}/initialallelefreqs.txt data");
$arg_hash->{allelefreqfile}="data/initialallelefreqs.txt";
$arg_hash->{initialhapmixlambdafile}="data/initiallambdas.txt";
$arg_hash->{fixedallelefreqs} = 0;
delete $arg_hash->{priorallelefreqfile};
#doAnalysis($executable,$arg_hash);
