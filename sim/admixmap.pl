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
    print "Results will be written to subdirectory $ENV{'RESULTSDIR'}";
    system($command);
    print "Starting R script to process output\n";
    system("RCMD BATCH --quiet --no-save --no-restore c:\\cvs\\genepi\\test\\AdmixmapOutput.R \
            $args->{resultsdir}/Rlog.txt");
    print "R script completed\n\n";
}

my $executable = 'c:/cvs/genepi/test/admixmap';

#########################################################################

my $arg_hash = 
{
    samples                    => 250, 
    burnin                     => 50,
    every                      => 5,
    numannealedruns            => 0,
    locusfile                  => 'data/loci.txt',
    genotypesfile              => 'data/genotypes.txt',
    priorallelefreqfile        => 'data/priorallelefreqs.txt',
    outcomevarfile             => 'data/outcome.txt',
    displaylevel               => 2,
# output files
    logfile                    => 'logfile.txt',
    paramfile                  => 'param.txt',
    regparamfile               => 'regparam.txt',
    hwscoretestfile            => 'HardyWeinbergTest.txt',
};

# model with prior allele freqs
$arg_hash->{resultsdir}            = 'priorAlleleFreqResults';  
$arg_hash->{indadmixturefile}   = "indivadmixture.txt";
$arg_hash->{dispersiontestfile}    = "dispersionTest.txt";
$arg_hash->{affectedsonlyscorefile} = 'affectedsonlyscores.txt';
$arg_hash->{ancestryassociationscorefile} = 'ancestryassociationscores.txt';
doAnalysis($executable,$arg_hash);
