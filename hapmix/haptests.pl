#!/usr/bin/perl
use strict; 
use File::Path;

################### DO NOT EDIT ABOVE THIS LINE ########################

# Change this to the location of the admixmap executable
my $executable = 'c:/cvs/genepi/test/admixmap';
my $arg_hash = 
{
#data files
    genotypesfile                   => 'data/genotypes.txt',
    locusfile                       => 'data/loci.txt',
    outcomevarfile                  => 'data/outcome.txt',
    
#main options
    samples  => 60,
    burnin   => 10,
    every    => 1,
    displaylevel => 3, 
    hapmixmodel => 1,
    
#output file options
    logfile                     => 'log.txt',
    paramfile                   => 'popadmixparams.txt',
    regparamfile                => 'regparam.txt',
    
# optional tests
    residualallelicassocscorefile => 'residualLDscoretests.txt',
    allelicassociationscorefile       => 'allelicassociationscoretests.txt',
    stratificationtestfile            => 'stratificationtest.txt',
    hwscoretestfile                   => 'HardyWeinbergtest.txt'
};

# model with reference prior on allele freqs in 2 populations
$arg_hash->{populations}           = 2;
$arg_hash->{resultsdir}            = 'TwoStateResults';  
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
    print "\nResults will be written to subdirectory $ENV{'RESULTSDIR'}";
    system($command);
    print "Starting R script to process output\n";
    system("RCMD BATCH --quiet --no-save --no-restore c:/cvs/genepi/test/AdmixmapOutput.R \
            $args->{resultsdir}/Rlog.txt");
    print "R script completed\n\n";
}



