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
my $executable = '../test/admixmap';
# command-line options are stored in an associative array (known as a hash in perl)  
my $arg_hash = {
#data files
    genotypesfile                   => 'data/genotypesAIMsOnly.txt',
    locusfile                       => 'data/lociAIMsOnly.txt',
#main options
    samples  => 1100,
    burnin   => 100,
    every    => 10,
    thermo   => 1,
    numannealedruns => 200,  
    displaylevel => 2,
#output file options
    logfile                     => 'log.txt',
};

# model with reference prior on allele freqs in 1 population, skin reflectance as continuous outcome var
$arg_hash->{populations}           = 1;
$arg_hash->{resultsdir}            = 'SinglePopResults';
#&doAnalysis($executable,$arg_hash);

# model with reference prior on allele freqs in 2 populations
$arg_hash->{populations}           = 2;
#$arg_hash->{samples}   = 6000;
#$arg_hash->{burnin}    = 1000;
$arg_hash->{paramfile}                 = 'popadmixparams.txt',
$arg_hash->{resultsdir}            = 'TwoPopsResults';  
&doAnalysis($executable,$arg_hash);

# model with reference prior on allele freqs in 3 populations
$arg_hash->{populations}           = 3;
$arg_hash->{resultsdir}            = 'ThreePopsResults';  
&doAnalysis($executable,$arg_hash);



