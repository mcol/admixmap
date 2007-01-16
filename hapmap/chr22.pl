#!/usr/bin/perl -w
# vim:set sw=4 softtabstop=4 ts=8 expandtab:
# 
# An example call (for a very fast run):
# ./chr22.pl --genotypes-file Eur/chr22data/genotypes5000_masked.txt \
# --locus-file Eur/chr22data/loci5000.txt --maskfile \
# Eur/chr22data/genotypes5000_index.txt --samples 12 --burnin 2 \
# --mutual-information

use strict;
use Getopt::Long;
use File::Path;
use Time::Local;

my $DEBUG = 0; # zero gives less output
print "script began: " . scalar(localtime()) . "\n";

##default values for script: Europeans, 4 states, 50 its with 10 burnin
my $samples=50;
my $burnin=10;
my $every=1;
my $POP = "Eur";
my $STATES = 4;
my $parallel = 0;
my $maskfile = '';
my $genotypes_file = '';
my $locus_file = '';
my $rerun = 0;
# Whether to calculate the mutual information
my $calculate_mi = 0;
my $usage = 0;

# Change this to the location of the hapmixmap executable
#my $executable = '../test/adm-para';
my $executable = '../test/hapmixmap';
#my $executable = '/ichec/home/users/doducd/test/hapmixmap';

##parse any command line options
GetOptions(
    "parallel!"  => \$parallel,
    "samples=i" => \$samples,
    "burnin=i"  => \$burnin,
    "every=i"   => \$every,
    "pop=s"     => \$POP,
    "states=i"  => \$STATES,
    "exec=s"    => \$executable,
    "mutual-information!"=> \$calculate_mi,
    "help!"     => \$usage,
    "genotypes-file=s"=> \$genotypes_file,
    "locus-file=s"  => \$locus_file,
    "re-run!"       => \$rerun,
    "maskfile=s"    => \$maskfile);

# print "Samples: $samples\n";
# exit(1);

# calculate_mutual_information();
# exit(0);

if ($usage) {
           ########################################################################
    print "\n";
    print "Usage: $0 [arguments]\n";
    print "  --parallel                Perform calculations in a parallel mode,\n";
    print "                            when on a cluster or multi-processor.\n";
    print "  --samples <integer>       Number of iterations\n";
    print "  --burnin <integer>        Lower than number of samples (above)\n";
    print "  --every <integer>         Integer, (samples - burnin) > (10 * every)\n";
    print "  --exec <file>             Name of the hapmixmap executable,\n";
    print "                            e.g. ../test/hapmixmap\n";
    print "  --pop [ Eur | Afr | Asi ] Population\n";
    print "  --states [integer]        Number of hidden states\n";
    print "  --mutual-information      Whether to calculate the mutual information\n";
    print "  --maskfile                File with indices of masked individuals\n";
    print "                            and loci\n";
    print "  --genotypes-file <name>   Genotypes data file\n";
    print "  --locus-file <name>       Locus data file\n";
    print "  --re-run                  Rerun with final values of lambda, freqs\n";
    print "                            in previous run as starting values\n";
    print "\n";
    exit(1);
}

my $datadir = "$POP/chr22data";
# $arg_hash is a hash of parameters passed to
# the executable as arguments.
#
# keys (left-hand side) are parameter names
# values (right-hand side) are parameter values
my $arg_hash = 
{
# data files
    # genotypesfile => "$datadir/genotypes5000_masked.txt",
    # locusfile     => "$datadir/loci5000.txt",

# phased data
    genotypesfile => "$datadir/genotypes_phased.txt",
    locusfile     => "$datadir/loci_phased.txt",

    #priorallelefreqfile             => 'data/priorallelefreqs.txt',
    #fixedallelefreqs => 1,
    populations=>$STATES,
    #outcomevarfile => 'chr22/dummyoutcome.txt',

checkdata=> 0,

#main options
    resultsdir => 'results',
    displaylevel   => 3, 

    samples  => $samples,
    burnin   => $burnin,
    every    => $every,

numannealedruns => 0,
thermo => 0,
hapmixmodel => 1,
#indadmixhiermodel => 0,
randommatingmodel => 0,

hapmixlambdaprior=>"30, 0.1, 10, 1",

allelefreqprior => "0.2, 1, 1",
#initialhapmixlambdafile => "$datadir/initialambdas.txt",
#allelefreqfile => "$datadir/initialallelefreqs.txt",

rhosamplerparams => "0.1, 0.00001, 10, 0.9, 20",

#output files
    logfile                     => 'logfile.txt',
    paramfile               => 'paramfile.txt',
    #regparamfile          => 'regparamfile.txt',
    #indadmixturefile     => 'indadmixture.txt',
    #ergodicaveragefile => 'ergodicaverage.txt',
    allelefreqprioroutputfile =>"allelefreqpriors.txt",
    allelefreqoutputfile  => "initialallelefreqs.txt",
    hapmixlambdaoutputfile => "initiallambdas.txt",

#optional tests
#residualallelicassocscorefile => 'residualLDscores.txt',
mhtestfile => "mhtest.txt",

    #allelicassociationscorefile       => 'allelicassociationscorefile.txt',
};

########################################################################
# If genotypes file and/or locus file was specified in the command-line,
# overwrite the arg_hash settings.
if ($genotypes_file) {
    $arg_hash->{'genotypesfile'} = $genotypes_file;
}
if ($locus_file) {
    $arg_hash->{'locusfile'} = $locus_file;
}

#model with $STATES block states

if (!$rerun) {
    # initial run
    $arg_hash->{resultsdir}="$POP/Results${STATES}States";
    doAnalysis($executable,$arg_hash);
} else {
    # rerun with final values of lambda, freqs in previous run as starting values
    $arg_hash->{resultsdir}="$POP/Results$STATES"."States";
    $arg_hash->{initialhapmixlambdafile} = "$datadir/initiallambdas.txt";
    $arg_hash->{allelefreqfile} = "$datadir/initialallelefreqs.txt";
    $arg_hash->{residualallelicassocscorefile} = 'residualLDscores.txt';
    doAnalysis($executable,$arg_hash);
}

#to gauge efficiency of Affymetrix chip
$arg_hash->{genotypesfile} = "$datadir/Affygenotypes.txt";
$arg_hash->{resultsdir} = "Affyresults";
#doAnalysis($executable,$arg_hash);

#to gauge efficiency of Illumina 300k chip
$arg_hash->{genotypesfile} = "$datadir/Ill300genotypes.txt";
$arg_hash->{resultsdir} = "Ill300results";
#doAnalysis($executable,$arg_hash);

#to gauge efficiency of Illumina 540k chip
$arg_hash->{genotypesfile} = "$datadir/Ill540genotypes.txt";
$arg_hash->{resultsdir} = "Ill540results";
#doAnalysis($executable,$arg_hash);

print "script ended: " . scalar(localtime()) . "\n";

sub getArguments
{
    my $hash = $_[0];
    my $filename = "args$POP$STATES.txt";
    open(OPTIONFILE, ">$filename") or die ("Could not open args file");
    foreach my $key (keys %$hash){
      print OPTIONFILE $key . '=' . $hash->{$key} . "\n";
    }

    # if we are using masked data
    if ($maskfile) {
        # If possible, append the contents of the index file with masked
        # loci to the option file, so users don't have to do it by hand.
        # Accept both full and $datadir-relative paths to the mask file.
        # $datadir-relative takes precedence.
        if (-r "$datadir/$maskfile") {
            $maskfile = "$datadir/$maskfile";
        }
        open(EXTERNAL_ARGS, "$maskfile")
            or warn("Can't open the external arguments file.");
        foreach my $line (<EXTERNAL_ARGS>) {
            print OPTIONFILE $line;
        }
        close(EXTERNAL_ARGS);
    }
    close OPTIONFILE;
    return " ".$filename;
}

sub calculate_mutual_information {
    my @r_call = qw(R CMD BATCH --no-save --no-restore);
    push(@r_call, "--population=$POP");
    push(@r_call, "--states=$STATES");
    push(@r_call, "MutualInformation.R");
    # "../test/AdmixmapOutput.R $args->{resultsdir}/Rlog.txt RESULTSDIR=$args->{resultsdir}");
    my $r_return = system(join(" ", @r_call));
    if ($r_return != 0) {
        print "R script returned an error code! $r_return\n";
        exit($r_return);
    }
}

sub runRscript
{
    my ($args) = @_;
    # Comment out the next three lines to run admixmap without R script
    print "Starting R script to process output\n";
    system("R CMD BATCH --quiet --no-save --no-restore ../test/AdmixmapOutput.R $args->{resultsdir}/Rlog.txt RESULTSDIR=$args->{resultsdir}");
    print "R script completed\n\n";
}

sub doAnalysis
{
    my ($prog,$args) = @_;
    my $command = "";
    if($parallel){
        $command = "mpiexec ";
        $args->{resultsdir} = $args->{resultsdir}. "Parallel";
    }
    $command = $command . $prog.getArguments($args);
    $ENV{'RESULTSDIR'} = $args->{resultsdir};
    if(system($command)==0){
        if($arg_hash->{allelefreqoutputfile}){
#copy file with final allele freqs to data dir ready for next time
            system("cp $arg_hash->{resultsdir}/$arg_hash->{allelefreqoutputfile} $datadir/$arg_hash->{allelefreqoutputfile}");
        }
        if($arg_hash->{hapmixlambdaoutputfile}){
#copy file with final lambdas to data dir ready for next time
            system("cp $arg_hash->{resultsdir}/$arg_hash->{hapmixlambdaoutputfile} $datadir/$arg_hash->{hapmixlambdaoutputfile}");
        }
        # Run the R script, change (1) to (0) to disable.
        if (1) {
            runRscript($args);
        }
        if ($calculate_mi) {
            calculate_mutual_information();
        }
    }
}
