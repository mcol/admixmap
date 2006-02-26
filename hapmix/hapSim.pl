#!/usr/bin/perl
use strict; 
use File::Path;

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

################### DO NOT EDIT ABOVE THIS LINE ########################
my $executable = 'c:/cvs/genepi/test/admixmap';

my $arg_hash = {
#data files
    genotypesfile                   => 'data/genotypes.txt',
    locusfile                       => 'data/loci.txt',
    #priorallelefreqfile             => 'data/allelefreqs.txt',
    #fixedallelefreqs                => 1,

    samples  => 250,
    burnin   => 50,
    every    => 1,
    numannealedruns => 0,
    displaylevel => 3, 
    globalsumintensitiesprior        => "40,1",
    logfile                          => 'log.txt',

    residualallelicassocscorefile    => 'residualLDTests.txt',
};

# model with single block states
$arg_hash->{populations}           = 1;
$arg_hash->{resultsdir}            = 'resultsSim1State';  
doAnalysis($executable,$arg_hash);

# model with 2 block states
$arg_hash->{hapmixmodel}           = 1;
$arg_hash->{populations}           = 2;
$arg_hash->{resultsdir}            = 'resultsSim2States';  
doAnalysis($executable,$arg_hash);

# model with 3 block states
$arg_hash->{populations}           = 3;
$arg_hash->{resultsdir}            = 'resultsSim3States';  
doAnalysis($executable,$arg_hash);

# model with 4 block states
$arg_hash->{populations}           = 4;
$arg_hash->{resultsdir}            = 'resultsSim4States';  
doAnalysis($executable,$arg_hash);

############### DO NOT EDIT BELOW THIS LINE ############################
