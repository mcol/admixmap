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
use File::Copy;
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
my $read_initial_params = 0;

# Change this to the location of the hapmixmap executable
#my $executable = '../test/adm-para';
my $executable = '/ichec/home/users/maciej/usr/bin/hapmixmap';
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
    print "                             Script will stop after masking.\n";
    print "                             Run it again with this option to do analysis.\n";
    print "  --percent-indivs <integer> Percent of masked individuals\n";
    print "  --percent-loci <integer>   Percent of masked loci\n";
    print "  --case-control-file <file> Genotypes file for case-control set\n";
    print "  --limit-loci <integer>     Limit number of loci used when masking.\n";
    print "\n";
    exit(1);
}

my %state_files = (
    arrivalratefile => "arrivalrates",
    mixturepropsfile => "mixtureprops",
    allelefreqfile => "allelefreqs",
    freqpriorfile => "freqprior",
);

my @archive_files = (
    'PPGenotypeProbs.txt',
    'EnergyTracePlot.ps',
    'logfile.txt',
    'loglikelihoodfile.txt',
    'PosteriorQuantiles.txt',
    'coefficient-of-constraint-by-locus-dput.txt',
    'coefficient-of-constraint-by-locus.txt',
    'mean-coefficient-of-constraint-no-uncert.txt',
    'mean-coefficient-of-constraint.txt',
);

my $datadir = "$POP/chr22data";
# $arg_hash is a hash of parameters passed to
# the executable as arguments.
#
# keys (left-hand side) are parameter names
# values (right-hand side) are parameter values
my $arg_hash = {
    deleteoldresults => 0,

#unphaseddata files
#    genotypesfile                   => "$datadir/genotypes5000.txt",
#    locusfile                       => "$datadir/loci5000.txt",

#phased data
    genotypesfile                   => "$datadir/phased_genotypes.txt",
    locusfile                       => "$datadir/phased_loci.txt",

#model
    states     => $STATES,
    checkdata       => 1,
    # mixturepropsprior => "50, 1",
    fixedmixtureprops => 1,
    fixedmixturepropsprecision =>1,

#main options
    resultsdir      => 'results',
    displaylevel    => 3,
    samples         => $samples,
    burnin          => $burnin,
    every           => $every,
    numannealedruns => 0,
    thermo          => 0,
    hapmixmodel     => 1,

#prior spec
    arrivalrateprior => "30, 0.1, 10, 1",
    # allelefreqprecisionprior => "0.2, 1, 1",
    allelefreqprecisionprior => "25, 1000, 10",

    arrivalratesamplerparams => "0.1, 0.00001, 10, 0.9, 20",

#output files
    logfile =>'logfile.txt',
    paramfile         =>'paramfile.txt',#mean and var of sampled arrival rates
    freqprecisionfile =>'freqprecision.txt', #mean and var of sampled allele freq precision
    arrivalrateposteriormeanfile => "ArrivalRatePosteriorMeans.txt",

#posterior means
    #arrivalrateposteriormeanfile => 'ArrivalRatePostMeans.txt',
#mean and var of sampled freq precision
    #allelefreqprecisionposteriormeanfile => 'FreqPrecisionPostMeans.txt',

    #ergodicaveragefile => 'ergodicaverage.txt',


#optional tests
    #residualallelicassocscorefile => 'residualLDscores.txt',
    # The following line breaks execution.
    # mhtestfile => "mhtest.txt",
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

# Data masking if requested. All the output files go to the ${datadir}
if ($mask_data) {
    print "Preparing train and test data.\n";

    my $prep_script = "cut-loci/mask-genotypes.py";
    my @data_prepare_args = (
        "python ${prep_script}",
        "--in-haplotypes ${datadir}/phased_genotypes.txt",
        "--in-loci ${datadir}/phased_loci.txt",
        "-i 17",
        "-l 30",
        "-m 5000",
        "--output-dir ${datadir}");

    my $dp_return = system(join(" ", @data_prepare_args));

    if ($dp_return != 0) {
        die("${prep_script} returned an error code.\n");
    }

    print "Masking finished.\n";
    print "Exiting after masking the data.\n";
    print "Please run the script without the --mask-data option\n";
    print "to perform the analysis.\n";
    exit(0);
}

########################################################################
# Analysis launch
########################################################################

$arg_hash->{resultsdir}="$POP/Chr22Results$STATES"."States2";

# Set the output state file names
my $opt_val;
foreach my $option_name (keys %state_files) {
    # Relative to the results directory.
    $arg_hash->{"final" . $option_name} = "state-" . $state_files{$option_name} . ".txt";
    if ($rerun) {
        # Relative to the working directory.
        $opt_val = $arg_hash->{resultsdir} . "/state-" . $state_files{$option_name} . "-latest.txt";
        # Give the initial file to the program only if the file exists.
        if (-r $opt_val) {
            $arg_hash->{"initial" . $option_name} = $opt_val;
        } else {
            warn("File $opt_val doesn't exist.\n");
        }
    }
}

doAnalysis($executable,$arg_hash);

print "script ended: " . scalar(localtime()) . "\n";

sub getArguments
{
    my $hash = $_[0];
    my $filename = "args$POP$STATES.txt";
    open(OPTIONFILE, ">$filename") or die ("Could not open args file");
    flock OPTIONFILE, 2;
    foreach my $key (sort keys %{$hash}){
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
    flock OPTIONFILE, 8;
    close OPTIONFILE;
    return $filename;
}

sub calculate_mutual_information {
    my @r_call = qw(R CMD BATCH --no-save --no-restore);
    push(@r_call, "--chromosome=Chr22");
    push(@r_call, "--population=$POP");
    push(@r_call, "--states=$STATES");
    push(@r_call, "MutualInformation.R");
    # "../test/AdmixmapOutput.R $args->{resultsdir}/Rlog.txt RESULTSDIR=$args->{resultsdir}");
    print "Calling the MutualInformation.R script.\n";
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
    system("R CMD BATCH --quiet --no-save --no-restore ../test/admixmap/AdmixmapOutput.R $args->{resultsdir}/Rlog.txt RESULTSDIR=$args->{resultsdir}");
    print "R script completed\n\n";
}

# Restore the saved train files. This should be done when interlacing
# training and testing runs.
sub restore_training_state {
    my $args = shift;
    my $slash = get_slash();
    for my $option_name (keys %state_files) {
        my $base_name = $args->{resultsdir} . $slash . "state-" . $state_files{$option_name};
        my $src_file = $base_name . "-trained-latest.txt";
        my $dst_file = $base_name . "-latest.txt";
        print "Copying $src_file to $dst_file.\n";
        copy($src_file, $dst_file)
            or die("Cannot copy '$src_file' to '$dst_file'.\n");
    }
}

sub get_slash {
    if($^O eq "MSWin32"){
        return "\\";
    } else {
        return "/";
    }
}

# After the program has finished, append a timestamp to all the
# output files and make a copy with "-latest" postfix for the
# next run.
sub rotate_files {
    my $args = shift;
    my $slash = get_slash();
    # Get the current timestamp in format YYYYMMDD-HHII
    my @rt = localtime();
    my $timestamp = sprintf("%04d%02d%02d-%02d%02d%02d", ($rt[5] + 1900), ($rt[4] + 1), $rt[3], $rt[2], $rt[1], $rt[0]);
    foreach my $option_name (keys %state_files) {
        my $base_name = $args->{resultsdir} . $slash . "state-" . $state_files{$option_name};
        my $src_file = $base_name . ".txt";
        my $timestamped = $base_name . "-" . $timestamp . ".txt";
        my $latest = $base_name . "-latest.txt";
        copy($src_file, $timestamped)
            or die("Cannot copy the '$src_file' to '$timestamped'.\n");
        # When training, archive the state.
        if (not $maskfile) {
            my @trained_files = (
                $base_name . "-trained-" . $timestamp . ".txt",
                $base_name . "-trained-latest.txt");
            foreach my $arc_file (@trained_files) {
                copy($src_file, $arc_file)
                    or die ("Cannot copy the '$src_file' to '$arc_file'.\n");
            }
        }
        move($src_file, $latest)
            or die("Cannot move the '$src_file' to '$latest'.\n");
    }
    foreach my $arc_file (@archive_files) {
        my $src_file = $args->{resultsdir} . "/" . $arc_file;
        # Chop off the extension
        my @parts = split(/\./, $src_file);
        my $ext = pop(@parts);
        my $dst_file = join(".", @parts) . "-$timestamp" . "." . $ext;
        if (-r $src_file) {
            copy($src_file, $dst_file) or die("Cannot copy '$src_file' to '$dst_file'.");
        } else {
            warn("Can't find '$src_file'.\n");
        }
    }
}

sub doAnalysis
{
    my ($prog,$args) = @_;
    my $command = "";
    if($parallel){
        $command = "mpiexec ";
        $args->{resultsdir} = $args->{resultsdir}. "parallel";
    }
    $command = $command . $prog . " " . getArguments($args);
    $ENV{'RESULTSDIR'} = $args->{resultsdir};
    if ($maskfile) {
        restore_training_state($args);
    }
    # Everybody likes Perl shortcuts.
    system($command) and die("'$command' has returned an error.\n");
    # Collect the time stats.
    system("bash log-time.sh " . $arg_hash->{resultsdir} . "/logfile.txt")
        and warn("Couldn't collect the time stats.");
    # run the r script
    runRscript($args);
    if ($calculate_mi) {
        calculate_mutual_information();
    }
    # Files should be rotated after mutual information calculation.
    rotate_files($args);
}
