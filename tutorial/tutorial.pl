#!/usr/bin/perl
use strict; 
use File::Path;

sub getArguments {
    my $hash = $_[0];
    my $filename = 'perlargs.txt';
    open(OPTIONFILE, ">$filename") or die ("Could not open args file");
    foreach my $key (keys %$hash){
      print OPTIONFILE $key . '=' . $hash->{$key} . "\n";
    }
    close OPTIONFILE;
    return " ".$filename;
}

sub doAnalysis {
    my ($prog,$args) = @_;
    my $command = $prog.getArguments($args);
    if (-e $args->{resultsdir}) {
	rmtree($args->{resultsdir});
    } 
    mkpath($args->{resultsdir});
    $ENV{'RESULTSDIR'} = $args->{resultsdir};
    print "\nResults will be written to subdirectory $ENV{'RESULTSDIR'}\n";
    system($command);
    my $rcmd = "R CMD";
    if($^O eq "MSWin32") {
	$rcmd = "Rcmd";
    }
    print "Starting R script to process output\n";
    system("$rcmd BATCH --quiet --no-save --no-restore ../test/AdmixmapOutput.R $args->{resultsdir}/Rlog.txt\n");
    print "R script completed\n\n";
}

################### DO NOT EDIT ABOVE THIS LINE ########################

# Change this to the location of the admixmap executable
my $executable = 'c:/cvs1/genepi/test/admixmap';
# command-line options are stored in an associative array (known as a hash in perl)  
my $arg_hash = {
#data files
    genotypesfile                   => 'data/genotypes.txt',
    locusfile                       => 'data/loci.txt',
    covariatesfile                  => 'data/covariates2std.txt', # age, sex 
    outcomevarfile                  => 'data/outcomevars.txt',
#main options
    samples  => 1200,
    burnin   => 200,
    every    => 5,
    numannealedruns => 0, #200, # 100, 
    displaylevel => 3,
#output file options
    logfile                     => 'log.txt',
    regparamfile                => 'regparam.txt',
    ergodicaveragefile          => 'cumulativeAverages.txt',
# optional tests
    #residualallelicassocscorefile       => 'residualLDscoretests.txt',
    #allelicassociationscorefile       => 'allelicassociationscoretests.txt',
    #haplotypeassociationscorefile     => 'hapassocscoretests.txt',
    stratificationtestfile            => 'stratificationtest.txt'
    #hwscoretestfile                   => 'HardyWeinbergTests.txt'
};

# model with reference prior on allele freqs in 1 population, skin reflectance as continuous outcome var
$arg_hash->{populations}           = 1;
$arg_hash->{resultsdir}            = 'SinglePopResults';
$arg_hash->{outcomes}              = 1,
$arg_hash->{targetindicator}       = 1; # skin reflectance
&doAnalysis($executable,$arg_hash);

# model with reference prior on allele freqs in 2 populations
$arg_hash->{populations}           = 2;
$arg_hash->{samples}   = 600;
$arg_hash->{burnin}    = 100;
$arg_hash->{every}    = 5;
$arg_hash->{paramfile}                 = 'popadmixparams.txt',
$arg_hash->{resultsdir}            = 'TwoPopsResults';  
&doAnalysis($executable,$arg_hash);

# model with reference prior on allele freqs in 3 populations
$arg_hash->{populations}           = 3;
$arg_hash->{resultsdir}            = 'ThreePopsResults';  
&doAnalysis($executable,$arg_hash);

# model with prior allele freqs 
delete $arg_hash->{populations};
$arg_hash->{resultsdir}                    = 'PriorFreqResultsSkin';  
$arg_hash->{samples}    = 1200;
$arg_hash->{burnin}    = 200;
$arg_hash->{priorallelefreqfile}           = 'data/priorallelefreqs.txt';
$arg_hash->{dispersiontestfile}            = 'dispersionTest.txt';
$arg_hash->{indadmixturefile}              = 'indivadmixture.txt';
$arg_hash->{ancestryassociationscorefile}  = 'ancestryassociationscorefile.txt';
&doAnalysis($executable,$arg_hash);

# model with prior allele freqs and diabetes as binary outcome var 
delete $arg_hash->{populations};
$arg_hash->{resultsdir}                = 'PriorFreqResultsDiabetes';  
$arg_hash->{targetindicator}           = 0; # diabetes as outcome
$arg_hash->{affectedsonlyscorefile}    = 'affectedsonlyscorefile.txt';
#$arg_hash->{thermo}    = 1;
#$arg_hash->{numannealedruns}    = 100;
&doAnalysis($executable,$arg_hash);

# model with historic allele freqs and both outcome vars
delete $arg_hash->{targetindicator};  
delete $arg_hash->{priorallelefreqfile};
delete $arg_hash->{dispersiontestfile};
delete $arg_hash->{affectedsonlyscorefile};
$arg_hash->{numannealedruns}    = 0;
$arg_hash->{resultsdir}                = 'HistoricFreqResults';  
$arg_hash->{randommatingmodel}         = 1;
$arg_hash->{outcomes}                  = 2; # both outcome vars - overrides targetindicator
$arg_hash->{historicallelefreqfile}    = "data/priorallelefreqs.txt";
$arg_hash->{etapriorfile}              = "data/etapriors.txt";
$arg_hash->{dispparamfile}             = "dispersionparams.txt";
$arg_hash->{fstoutputfile}             = "lociFst.txt";
$arg_hash->{allelefreqoutputfile}      = "allelefreqs.txt";
&doAnalysis($executable,$arg_hash);

