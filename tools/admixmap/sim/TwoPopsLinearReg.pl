#!/usr/bin/perl
use strict; 
use File::Path;

sub getArguments
{
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

    $ENV{'RESULTSDIR'} = $args->{resultsdir};
    print "\nResults will be written to subdirectory $ENV{'RESULTSDIR'}\n";
    system($command);
    my $rcmd = "R CMD";
    if($^O eq "MSWin32") {
	$rcmd = "Rcmd";
    }
    print "Starting R script to process output\n";
    system("copy LocusTable.txt $args->{resultsdir}\n");
    system("$rcmd BATCH --quiet --no-save --no-restore ..\\..\\..\\dist\\admixmap\\AdmixmapOutput.R $args->{resultsdir}/Rlog.txt\n");
    print "R script completed\n\n";
}

################### DO NOT EDIT ABOVE THIS LINE ########################

# Change this to the location of the admixmap executable
my $executable = '../../src/admixmap/admixmap';
if($^O eq "MSWin32") {
    $executable = '..\\..\\..\\src\\admixmap\\admixmap';
}
#$executable = "s:\\sharedfolders\\genepi\\admix\\admixmapbackups\\archive\\admixmap3.2b";
$executable = "c:\\oldgenepi\\trunk\\test\\admixmap";

# command-line options are stored in an associative array (known as a hash in perl)  
my $arg_hash = {
#data files
    genotypesfile                   => 'data/genotypes.txt',
    locusfile                       => 'data/loci.txt',
    outcomevarfile                  => 'data/outcome.txt',
#main options
    samples  => 6000,
    burnin   => 1000,
    every    => 10,
    numannealedruns => 0, #200, # 100, 
    displaylevel => 2,
#output file options
    logfile                     => 'log.txt',
    regparamfile                => 'regparam.txt',
};

# model with reference prior on allele freqs in 2 populations
$arg_hash->{populations}           = 2;
$arg_hash->{paramfile}                 = 'popadmixparams.txt',
$arg_hash->{resultsdir}            = 'TwoPopsResults';  
&doAnalysis($executable,$arg_hash);



