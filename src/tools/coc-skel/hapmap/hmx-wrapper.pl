#!/usr/bin/perl -w
# vim:set sw=4 softtabstop=4 ts=8 expandtab:
# 
# An example call (for a very fast run):
# ./hmx-wrapper.pl --genotypes-file Eur/chr22data/genotypes5000_masked.txt \
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
use File::Temp qw/ :POSIX /;

my $DEBUG = 0; # zero gives less output

##default values for script: Europeans, 4 states, 50 its with 10 burnin
my $samples=12;
my $burnin=1;
my $every=1;
my $POP = "";
my $STATES = 8;
my $genotypes_file = '';
my $locus_file = '';
my $rerun = 0;
# Whether to calculate the mutual information
my $usage = 0;
# Whether to perform data masking first, for the phased data
my $read_initial_params = 0;
my $chromosome = "";

# Change this to the location of the hapmixmap executable
#my $executable = '../test/adm-para';
my $executable = $ENV{'HOME'} . "/usr/bin/hapmixmap";
#my $executable = '/ichec/home/users/doducd/test/hapmixmap';

##parse any command line options
GetOptions(
    "samples=i" => \$samples,
    "chromosome=i" => \$chromosome,
    "burnin=i"  => \$burnin,
    "every=i"   => \$every,
    "pop=s"     => \$POP,
    "states=i"  => \$STATES,
    "exec=s"    => \$executable,
    "help!"     => \$usage,
    "genotypes-file=s"  => \$genotypes_file,
    "locus-file=s"      => \$locus_file,
    "re-run!"           => \$rerun);

if (not $chromosome) {
    print "Please specify the chromosome.\n";
    $usage = 1;
}

if (not $POP) {
    print "Please specify population.\n";
    $usage = 1;
}

if ($usage) {
    print "\n";
    print "Usage: $0 [arguments]\n";
    print "  --samples <integer>        Number of iterations\n";
    print "  --burnin <integer>         Lower than number of samples (above)\n";
    print "  --every <integer>          Integer, (samples - burnin) > (10 * every)\n";
    print "  --exec <file>              Name of the hapmixmap executable,\n";
    print "                             e.g. ../test/hapmixmap\n";
    print "  --pop [ ceu | yri | jpt-chb ]  Population\n";
    print "  --states <integer>         Number of hidden states\n";
    print "  --genotypes-file <name>    Genotypes data file\n";
    print "  --locus-file <name>        Locus data file\n";
    print "  --re-run                   Rerun with final values of lambda, freqs\n";
    print "                             in previous run as starting values\n";
    print "\n";
    exit(1);
}

# Options are checked, we're free to go.
print "script began: " . scalar(localtime()) . "\n";

# The final state file name options are absent in hapmixmap, the
# following state files saving doesn't work anymore.
my %state_files = (
    # arrivalratefile => "arrivalrates",
    # mixturepropsfile => "mixtureprops",
    # allelefreqfile => "allelefreqs",
    # freqpriorfile => "freqprior",
);

my @archive_files = (
    'PPGenotypeProbs.txt',
    'EnergyTracePlot.ps',
    'logfile.txt',
    'args.txt',
    'loglikelihoodfile.txt',
    'PosteriorQuantiles.txt',
    'coefficient-of-constraint-by-locus-dput.txt',
    'coefficient-of-constraint-by-locus.txt',
    'mean-coefficient-of-constraint-no-uncert.txt',
    'mean-coefficient-of-constraint.txt',
);

my $datadir = "hapmap-hapmixmap";
my $padded_chr = sprintf("%02d", $chromosome);

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
    # genotypesfile                   => "$datadir/chr${chromosome}-${POP}-genotypes.txt",
    # locusfile                       => "$datadir/chr${chromosome}-${POP}-loci.txt",

    # FPHD-formatted data
    genotypesfile                   => "$datadir/${POP}-chr${padded_chr}-gt.txt",
    locusfile                       => "$datadir/${POP}-chr${padded_chr}-loci.txt",

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
    arrivalrateprior => "1.2, 0.05, 1, 0.5",
    allelefreqprecisionprior => "0.25, 1",
    freqprecisionhiermodel => 0,

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

########################################################################
# Analysis launch
########################################################################

$arg_hash->{resultsdir}="results/chr${chromosome}-${POP}";

# # Set the output state file names
# my $opt_val;
# foreach my $option_name (keys %state_files) {
#     # Relative to the results directory.
#     $arg_hash->{"final" . $option_name} = "state-" . $state_files{$option_name} . ".txt";
#     if ($rerun) {
#         # Relative to the working directory.
#         $opt_val = $arg_hash->{resultsdir} . "/state-" . $state_files{$option_name} . "-latest.txt";
#         # Give the initial file to the program only if the file exists.
#         if (-r $opt_val) {
#             $arg_hash->{"initial" . $option_name} = $opt_val;
#         } else {
#             warn("File $opt_val doesn't exist.\n");
#         }
#     }
# }

if ($rerun) {
    $arg_hash->{"initialvaluedir"} = $arg_hash->{resultsdir};
}

doAnalysis($executable,$arg_hash);

# #to gauge efficiency of Affymetrix chip
# $arg_hash->{genotypesfile} = "$datadir/Affygenotypes.txt";
# $arg_hash->{resultsdir} = "Affyresults";
# #doAnalysis($executable,$arg_hash);
# 
# #to gauge efficiency of Illumina 300k chip
# $arg_hash->{genotypesfile} = "$datadir/Ill300genotypes.txt";
# $arg_hash->{resultsdir} = "Ill300results";
# #doAnalysis($executable,$arg_hash);
# 
# #to gauge efficiency of Illumina 540k chip
# $arg_hash->{genotypesfile} = "$datadir/Ill540genotypes.txt";
# $arg_hash->{resultsdir} = "Ill540results";
# #doAnalysis($executable,$arg_hash);

print "script ended: " . scalar(localtime()) . "\n";

sub getArguments
{
    my $hash = $_[0];
    my $filename = tmpnam();
    open(OPTIONFILE, ">$filename") or die ("Could not open args file");
    foreach my $key (sort keys %{$hash}){
      print OPTIONFILE $key . '=' . $hash->{$key} . "\n";
    }

    close OPTIONFILE;
    return $filename;
}

sub calculate_mutual_information {
    my @r_call = qw(R CMD BATCH --no-save --no-restore);
    push(@r_call, "--chromosome=Chr$chromosome");
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
        if (1) {
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
    my $extra_options = '';
    my $extra_options_file = "extra-cli-options.txt";
    if (-r $extra_options_file) {
        open(EXTRA_OPTS, "<$extra_options_file");
        my @tmp = <EXTRA_OPTS>;
        $extra_options = $tmp[0];
        close(EXTRA_OPTS);
    }
    my $option_pfx = '';
    if ($extra_options) {
        $option_pfx = '-f';
    }
    my $args_file_name = getArguments($args);
    $command = $command . $prog . " " . $option_pfx . $args_file_name . " $extra_options";
    $ENV{'RESULTSDIR'} = $args->{resultsdir};
    # Everybody likes Perl shortcuts.
    system($command) and die("'$command' has returned an error.\n");
    # Collect the time stats.
    system("bash log-time.sh >> time-stats.txt " . $arg_hash->{resultsdir} . "/logfile.txt")
        and warn("Couldn't collect the time stats.");
    # run the r script
    runRscript($args);
    # Files should be rotated after mutual information calculation.
    rotate_files($args);
    unlink($args_file_name); # clean up temporary file
}

