#!/usr/bin/perl
my $doreturn = do "../../../dist/doanalysis.pl";
unless($doreturn) {
    die "couldn't parse doanalysis.pl: $@" if $@;
    #die "couldn't read doanalysis.pl: $!"    unless defined $doreturn;
    #die "couldn't run doanalysis.pl"       unless $doreturn;
}

################### DO NOT EDIT ABOVE THIS LINE ########################

# Change these to the locations of the admixmap executable and R script
my $prog = '../../../src/admixmap/admixmap';
my $rscript = "../../../dist/admixmap/AdmixmapOutput.R";
#$prog = "s:\\sharedfolders\\genepi\\admix\\admixmapbackups\\archive\\admixmap3.2b";
#$prog = "c:\\oldgenepi\\trunk\\test\\admixmap";

my $arg_hash = { # command-line options are stored in a hash (associative array)  
#data files
    genotypesfile                   => 'data/genotypes.txt',
    locusfile                       => 'data/loci.txt',
    outcomevarfile                  => 'data/outcome.txt',
#main options
    samples  => 1100,
    burnin   => 100,
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
$arg_hash->{regressionpriorprecision}            = 0.01;
$arg_hash->{priorallelefreqfile}            = "data/trueallelefreqs.txt";
$arg_hash->{fixedallelefreqs}            = 1;
doAnalysis($prog, $rscript, $arg_hash);
copy("$arg_hash->{resultsdir}/PosteriorQuantiles.txt", "PosteriorQuantiles1998.txt") 
    or die "File cannot be copied";

$prog = "s:\\sharedfolders\\genepi\\admix\\admixmapbackups\\archive\\admixmap3.2b";
delete $arg_hash->{regressionpriorprecision};
#&doAnalysis($prog, $rscript, $arg_hash);
#copy("$arg_hash->{resultsdir}/PosteriorQuantiles.txt", "PosteriorQuantiles3.2b.txt")
#    or die "File cannot be copied";



