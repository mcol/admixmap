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
    system("R CMD BATCH --quiet --no-save --no-restore ~/genepi/trunk/dist/AdmixmapOutput.R $args->{resultsdir}/Rlog.txt");
    print "R script completed\n\n";
}

my $executable = 'admixmap';

#########################################################################

my $arg_hash = 
{
    samples                    => 650, 
    burnin                     => 150,
    every                      => 5,
    thermo                     => 0,
    numannealedruns            => 0,
    locusfile                  => "data/loci.txt",
    genotypesfile              => "data/genotypes.txt",
    outcomevarfile             => 'data/outcome.txt',
    #populations                => 2,
    displaylevel              => 3,
    
# output files
    logfile                    => 'logfile.txt',
    allelicassociationscorefile => 'allelicassociationscores.txt'
};

# model with prior allele freqs
$arg_hash->{fixedallelefreqs}      = 0;
$arg_hash->{priorallelefreqfile}   = "data/priorallelefreqs.txt";
$arg_hash->{resultsdir}            = 'results';  
$arg_hash->{indadmixturefile}   = "indivadmixture.txt";
$arg_hash->{dispersiontest}    = 1;
$arg_hash->{affectedsonlytest} = 1;
$arg_hash->{ancestryassociationtest} = 1;
$arg_hash->{hwtest} = 1;
$arg_hash->{residualldtest} = 1;
doAnalysis($executable,$arg_hash);

