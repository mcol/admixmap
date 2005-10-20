#!/usr/bin/perl
use strict; 
use File::Path;

################### DO NOT EDIT ABOVE THIS LINE ########################

# Change this to the location of the admixmap executable
my $executable = '../test/admixmap.exe';
# command-line options are stored in an associative array (known as a hash in perl)  
my $arg_hash = 
{
#data files
    genotypesfile                   => 'data/genotypes.txt',
    locusfile                       => 'data/loci.txt',
    #covariatesfile                  => 'data/covariates3std.txt', # age, sex, income
    covariatesfile                  => 'data/covariates2std.txt', # age, sex only
    outcomevarfile                  => 'data/outcomevars.txt',
    
#main options
    samples  => 1200,
    burnin   => 200,
    every    => 5,
    coutindicator => 0, 
    
#output file options
    logfile                     => 'log.txt',
    paramfile                   => 'popadmixparams.txt',
    regparamfile                => 'regparam.txt',
    
# optional tests
    allelicassociationscorefile       => 'allelicassociationscoretests.txt',
    haplotypeassociationscorefile     => 'hapassocscoretests.txt',
    stratificationtestfile            => 'stratificationtest.txt',
    hwscoretestfile                   => 'HardyWeinbergTest.txt'
};

# model with reference prior on allele freqs in 1 population, skin reflectance as outcome var
$arg_hash->{populations}           = 1;
$arg_hash->{resultsdir}            = 'SinglePopResults';  
$arg_hash->{analysistypeindicator} = 2; # continuous outcome var
$arg_hash->{targetindicator}       = 1; # skin reflectance
doAnalysis($executable,$arg_hash);

# model with reference prior on allele freqs in 2 populations
$arg_hash->{populations}           = 2;
$arg_hash->{resultsdir}            = 'TwoPopsResults';  
#$arg_hash->{samples}               = 5500;
#$arg_hash->{burnin}                = 500;
#$arg_hash->{every}                 = 5;
doAnalysis($executable,$arg_hash);

# model with reference prior on allele freqs in 3 populations
#$arg_hash->{populations}           = 3;
$arg_hash->{resultsdir}            = 'ThreePopsResults';  
doAnalysis($executable,$arg_hash);

# model with prior allele freqs 
delete $arg_hash->{populations};
$arg_hash->{resultsdir}                    = 'PriorFreqResultsSkin';  
$arg_hash->{priorallelefreqfile}           = 'data/priorallelefreqs.txt';
$arg_hash->{dispersiontestfile}            = 'dispersionTest.txt';
$arg_hash->{indadmixturefile}              = 'indivadmixture.txt';
$arg_hash->{ancestryassociationscorefile}  = 'ancestryassociationscorefile.txt';
#$arg_hash->{samples}               = 1100;
#$arg_hash->{burnin}                = 100;
#$arg_hash->{every}                 = 5;
doAnalysis($executable,$arg_hash);

# model with prior allele freqs and diabetes as outcome var 
delete $arg_hash->{populations};
$arg_hash->{resultsdir}                = 'PriorFreqResultsDiabetes';  
$arg_hash->{analysistypeindicator}     = 3; # binary outcome var
$arg_hash->{targetindicator}           = 0; # diabetes
$arg_hash->{affectedsonlyscorefile}    = 'affectedsonlyscorefile.txt';
doAnalysis($executable,$arg_hash);

# model with historic allele freqs and both outcome vars
delete $arg_hash->{priorallelefreqfile};
delete $arg_hash->{dispersiontestfile};
delete $arg_hash->{affectedsonlyscorefile};
$arg_hash->{resultsdir}                = 'HistoricFreqResults';  
$arg_hash->{randommatingmodel}         = 1;
$arg_hash->{analysistypeindicator}     = 5; # both outcome vars
$arg_hash->{historicallelefreqfile}    = "data/priorallelefreqs.txt";
$arg_hash->{etapriorfile}              = "data/etapriors.txt";
$arg_hash->{dispparamfile}             = "dispersionparams.txt";
$arg_hash->{fstoutputfile}             = "lociFst.txt";
$arg_hash->{allelefreqoutputfile}      = "allelefreqs.txt";
doAnalysis($executable,$arg_hash);

############### DO NOT EDIT BELOW THIS LINE ############################

sub getArguments
{
    my $hash = $_[0];
    my $arg = '';
    foreach my $key (keys %$hash){
	$arg .= ' --'. $key .'='. $hash->{$key};
    }
    return $arg;
}

sub doAnalysis
{
    my ($prog,$args) = @_;
    my $command = $prog.getArguments($args);
    if (-e $args->{resultsdir}) {
	rmtree($args->{resultsdir});
    } 
    mkpath($args->{resultsdir});
    $ENV{'RESULTSDIR'} = $args->{resultsdir};
    print "Results will be written to subdirectory $ENV{'RESULTSDIR'}\n";
    system($command);
    print "Starting R script to process output\n";
    system("RCMD BATCH --quiet --no-save --no-restore ../test/AdmixmapOutput.R $args->{resultsdir}/Rlog.txt");
    print "R script completed\n\n";
}




