#!/usr/bin/perl
use strict; 
use File::Path;
use Getopt::Std;



#---------------------------------------------------------------------------
# Command-line options and usage errors:
#---------------------------------------------------------------------------

my $PROG = ($0 =~ m/.*?([^\/]+)$/)[ 0 ];

sub usage
    {
    print STDERR "\nUsage: $PROG [-h] [-p] [-x <admixmap-executable>]\n\n";
    exit 1
    }


my %args;
getopts( "x:hp", \%args ) or usage;

usage if defined $args{"h"};

my $executable = defined $args{'x'} ? $args{'x'} : 'admixmap';
my $usepedfile = defined $args{'p'};



#---------------------------------------------------------------------------
# getArguments
#---------------------------------------------------------------------------

sub getArguments
{
    my $hash = $_[0];
    my $arg = '';
    foreach my $key (keys %$hash){
	$arg .= ' --'. $key .'='. $hash->{$key};
    }
    return $arg;
}



#---------------------------------------------------------------------------
# doAnalysis
#---------------------------------------------------------------------------

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
    print "Command: " . $command . "\n";
    my $status = system($command);
    if ( $status == -1 )
	{
	print STDERR $PROG . ": error executing: $command : $!\n";
	}
    elsif ( $status != 0 )
	{
	print STDERR $PROG . ": $prog exited with non-zero status: $status\n";
	}

    print "Starting R script to process output\n";
    system("R CMD BATCH --quiet --no-save --no-restore ~/genepi/trunk/dist/AdmixmapOutput.R $args->{resultsdir}/Rlog.txt");
    print "R script completed\n\n";
}



#########################################################################

my $arg_hash = 
{
    samples                    => 650, 
    burnin                     => 150,
    every                      => 5,
    thermo                     => 0,
    numannealedruns            => 0,
    locusfile                  => "data/loci.txt",
    genotypesfile              => ($usepedfile ? "data/genotypes.ped" : "data/genotypes.txt"),
    outcomevarfile             => 'data/outcome.txt',
    #populations                => 2,
    displaylevel              => 3,
    
# output files
    logfile                    => 'logfile.txt',
    allelicassociationscorefile => 'allelicassociationscores.txt'
};


# **** These should be re-enabled when the are supported by pedigrees: ****
if ( $usepedfile ) {
    $arg_hash->{"no-conjugate-update"} = 1;
    $arg_hash->{fixedallelefreqs} = 1;
    $arg_hash->{residualldtest} = 0;
} else {
    $arg_hash->{fixedallelefreqs} = 0;
    $arg_hash->{residualldtest} = 1;
}


# model with prior allele freqs
$arg_hash->{priorallelefreqfile}   = "data/priorallelefreqs.txt";
$arg_hash->{resultsdir}            = 'results';  
$arg_hash->{indadmixturefile}   = "indivadmixture.txt";
$arg_hash->{dispersiontest}    = 1;
$arg_hash->{affectedsonlytest} = 1;
$arg_hash->{ancestryassociationtest} = 1;
$arg_hash->{hwtest} = 1;
doAnalysis($executable,$arg_hash);

