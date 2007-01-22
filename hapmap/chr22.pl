#!/usr/bin/perl -w
# vim:set sw=4 softtabstop=4 ts=8 expandtab:
# 
# An example call (for a very fast run):
# ./chr22.pl --genotypes-file Eur/chr22data/genotypes5000_masked.txt \
# --locus-file Eur/chr22data/loci5000.txt \
# --maskfile Eur/chr22data/genotypes5000_index.txt \
# --samples 12 --burnin 2 \
# --mutual-information
#
# TODO: Filenames should be standarized. There are still some hard-coded
# file names which may be confusing to users.
#
# Mutual information calculations. I don't know if using a Perl script
# as a notebook is great, but I'll do it anyway.
#
# Source files:
# phased_genotypes.txt (haploid)
# phased_loci.txt
#
# Those files are going to be split into three files. Number of loci
# will be trimmed to 5000. File names:
# mi_genotypes.txt    (haploid) (50 individuals)
# mi_cc_genotypes.txt (diploid) (10 individuals)
# mi_loci.txt
#
# The file mi_cc_genotypes.txt is going to be masked, all
# individuals, in R.

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
# Whether to perform data masking first, for the phased data
my $mask_data = 0;
my $mask_percent_indivs = 0;
my $mask_percent_loci = 0;
my $ccgenotypesfile = '';
my $limit_loci = 0;

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
    "maskfile=s"    => \$maskfile,
    "mask-data!"     => \$mask_data,
    "percent-indivs=i" => \$mask_percent_indivs,
    "percent-loci=i"  => \$mask_percent_loci,
    "case-control-file=s" => \$ccgenotypesfile,
    "limit-loci=i" => \$limit_loci);

if ($mask_data) {
    if ($mask_percent_indivs <= 0 || $mask_percent_indivs > 100) {
        print "Percent of masked individuals must be between 1 and 100.\n";
        exit(1);
    }
    if ($mask_percent_loci <= 0 || $mask_percent_loci > 100) {
        print "Percent of masked loci must be between 1 and 100.\n";
        exit(1);
    }
}

if ($usage) {
           ########################################################################
    print "\n";
    print "Usage: $0 [arguments]\n";
    print "  --parallel                 Perform calculations in a parallel mode,\n";
    print "                             when on a cluster or multi-processor.\n";
    print "  --samples <integer>        Number of iterations\n";
    print "  --burnin <integer>         Lower than number of samples (above)\n";
    print "  --every <integer>          Integer, (samples - burnin) > (10 * every)\n";
    print "  --exec <file>              Name of the hapmixmap executable,\n";
    print "                             e.g. ../test/hapmixmap\n";
    print "  --pop [ Eur | Afr | Asian ]  Population\n";
    print "  --states [integer]         Number of hidden states\n";
    print "  --mutual-information       Whether to calculate the mutual information\n";
    print "  --maskfile                 File with indices of masked individuals\n";
    print "                             and loci\n";
    print "  --genotypes-file <name>    Genotypes data file\n";
    print "  --locus-file <name>        Locus data file\n";
    print "  --re-run                   Rerun with final values of lambda, freqs\n";
    print "                             in previous run as starting values\n";
    print "  --mask-data                Whether to perform data masking first.\n";
    print "                             In such a case, genotypes_phased.txt\n";
    print "                             and loci_phased.txt will be processed.\n";
    print "  --percent-indivs <integer> Percent of masked individuals\n";
    print "  --percent-loci <integer>   Percent of masked loci\n";
    print "  --case-control-file <file> Genotypes file for case-control set\n";
    print "  --limit-loci <integer>     Limit number of loci used when masking.\n";
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

#unphaseddata files
#    genotypesfile                   => "$datadir/genotypes5000.txt",
#    locusfile                       => "$datadir/loci5000.txt",

#phased data
    genotypesfile                   => "$datadir/phased_genotypes.txt",
    locusfile                       => "$datadir/phased_loci.txt",

    #priorallelefreqfile => 'data/priorallelefreqs.txt',
    #fixedallelefreqs => 1,
    states     => $STATES,
    #outcomevarfile => 'chr22/dummyoutcome.txt',
    checkdata       => 0,

#main options
    resultsdir      => 'results',
    displaylevel    => 3,
    samples         => $samples,
    burnin          => $burnin,
    every           => $every,
    numannealedruns => 0,
    thermo          => 0,
    hapmixmodel     => 1,
    randommatingmodel => 0,

#prior spec
    hapmixlambdaprior => "30, 0.1, 10, 1",
    allelefreqprior => "0.2, 1, 1",

#initial values
    #initialhapmixlambdafile => "$datadir/initialambdas.txt",
    #allelefreqfile => "$datadir/initialallelefreqs.txt",
    rhosamplerparams => "0.1, 0.00001, 10, 0.9, 20",
#output files
    logfile =>'logfile.txt',
    paramfile =>'paramfile.txt',#mean and var of sampled arrival rates
    dispparamfile => "allelefreqpriors.txt",#mean and var of sampled freq dispersion

    #regparamfile          => 'regparamfile.txt',
    #ergodicaveragefile => 'ergodicaverage.txt',

    allelefreqprioroutputfile   => "initialfreqdispersion.txt",
    allelefreqoutputfile        => "initialallelefreqs.txt",
    hapmixlambdaoutputfile      => "initiallambdas.txt",

#optional tests
    #residualallelicassocscorefile => 'residualLDscores.txt',
    mhtestfile => "mhtest.txt",
    #allelicassociationscorefile       => 'allelicassociationscorefile.txt',
};

# If genotypes file and/or locus file was specified in the command-line,
# overwrite the arg_hash settings.
if ($genotypes_file) {
    $arg_hash->{'genotypesfile'} = "$datadir/$genotypes_file";
}
if ($locus_file) {
    $arg_hash->{'locusfile'} = "$datadir/$locus_file";
}
if ($ccgenotypesfile) {
    $arg_hash->{'ccgenotypesfile'} = "$datadir/$ccgenotypesfile";
}

#model with $STATES block states


# Data masking if requested
# Data set will be base-named "train", that is:
# masked_genotypes.txt
# masked_loci.txt
#
#
if ($mask_data) {
    my $train_basename = "mi";
    print "Preparing train and test data.\n";
    my $lf = $arg_hash->{'locusfile'};
    my $new_locus_file = "mi_loci.txt";
    my @data_prepare_args = (
        "./prepareTestingData.pl ",
        "--phased-file $datadir/phased_genotypes.txt ",
        "--haploid-file $datadir/${train_basename}_genotypes.txt ",
        "--case-control-file $datadir/${train_basename}_cc.txt ",
        "--percent-indivs $mask_percent_indivs ",
        "--limit-loci $limit_loci ",
        "--in-locus-file $lf ",
        "--out-locus-file $datadir/$new_locus_file");
    my $dp_return = system(join(" ", @data_prepare_args));
    $arg_hash->{'locusfile'} = "$datadir/$new_locus_file";
    if ($dp_return != 0) {
        die("prepareTestingData.pl returned an error code.\n");
    }
    # Count the lines in the train data to determine the number of
    # individuals
    my $train_file_name = "$datadir/${train_basename}_genotypes.txt";
    open(TRAIN_DATA, "<$train_file_name")
        or die("Couldn't open the train data file: '$train_file_name'.");
    # One line for headers and two lines per individual
    my @train_array = <TRAIN_DATA>;
    my $train_indivs = (scalar(@train_array) - 1) / 2;
    close TRAIN_DATA;
    print "Masking data: $mask_percent_indivs% of individuals";
    print " and $mask_percent_loci% of loci. Running R for this.\n";
    print "Train individuals: $train_indivs.\n";
    my @r_call = qw(R CMD BATCH --no-save --no-restore);
    push(@r_call, "--population=$POP");
    push(@r_call, "--basename=${train_basename}_cc");
    push(@r_call, "--loci-file=$datadir/$locus_file");
    # All the individuals in the case-control file should be masked,
    # hence 100% masked individuals.
    push(@r_call, "--percent-indivs=100");
    push(@r_call, "--percent-loci=$mask_percent_loci");
    push(@r_call, "--indiv-offset=$train_indivs");
    push(@r_call, "maskGenotypes.R");
    my $r_return = system(join(" ", @r_call));
    if ($r_return != 0) {
        die "maskGenotypes.R script returned an error code: $r_return\n";
    }
    print "Masking finished.\n";
    print "Exiting after masking the data.\n";
    print "Please run the script without --mask-data option\n";
    print "to perform the analysis.\n";
    exit(0);
}

########################################################################
# Analysis launch
########################################################################

if (!$rerun) {
    # initial run
    $arg_hash->{resultsdir}="$POP/Chr22Results$STATES"."States2";
    doAnalysis($executable,$arg_hash);
} else {
    # rerun with final values of lambda, freqs, freq dispersion in previous run as starting values
    $arg_hash->{resultsdir}="$POP/Chr22Results$STATES"."States2";
    $arg_hash->{initialhapmixlambdafile} = "$datadir/initiallambdas.txt";
    $arg_hash->{allelefreqfile} = "$datadir/initialallelefreqs.txt";
    $arg_hash->{initialfreqdispersionfile} = "$datadir/initialfreqdispersion.txt";
    #$arg_hash->{allelicassociationscorefile} = 'AllelicAssocScores.txt';
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
    push(@r_call, "--chromosome=Chr22");
    push(@r_call, "--population=$POP");
    push(@r_call, "--states=$STATES");
    push(@r_call, "MutualInformation.R");
    # "../test/AdmixmapOutput.R $args->{resultsdir}/Rlog.txt RESULTSDIR=$args->{resultsdir}");
    my $r_return = system(join(" ", @r_call));
    if ($r_return != 0) {
        die("MutualInformation.R script returned an error code! $r_return\n");
    }
}

# Run R script which processes output files.
sub runRscript
{
    my ($args) = @_;
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
        $args->{resultsdir} = $args->{resultsdir}. "parallel";
    }
    $command = $command . $prog.getArguments($args);
    $ENV{'RESULTSDIR'} = $args->{resultsdir};
    my $returncode = system($command);
    if($returncode == 0){
##program ran successfully
#copy initial value files from resultsdir to datadir
#using wildcard to copy all in one go
        my $copycmd = "cp";
        my $slash = "/";
        if($^O eq "MSWin32"){
            $copycmd = "copy";
            $slash = "\\";
        }
        system("$copycmd $args->{resultsdir}$slash"."initial*.txt $datadir");
## run the r script
        runRscript($args);
        if ($calculate_mi) {
            calculate_mutual_information();
        }
    }
    else{
        warn("hapmixmap has returned an error code $returncode.");
    }
}
