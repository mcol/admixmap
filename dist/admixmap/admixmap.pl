#!/usr/bin/perl -w
my $doreturn = do "../doanalysis.pl";
unless($doreturn) {
    die "couldn't parse doanalysis.pl: $@" if $@;
    #die "couldn't read doanalysis.pl: $!"    unless defined $doreturn;
    #die "couldn't run doanalysis.pl"       unless $doreturn;
}

# Change these filepaths to the location of the program and R script
my $prog = '../../src/admixmap/admixmap';
my $rscript = './AdmixmapOutput.R';

# $arg_hash is a hash of parameters passed to
# the executable as arguments.
#
# keys (left-hand side) are parameter names
# values (right-hand side) are parameter values
my $arg_hash = {
 ##data files
 genotypesfile        => 'data/genotypes.txt',
 locusfile            => 'data/loci.txt',
 priorallelefreqfile  => 'data/priorallelefreqs.txt',
 #populations => 1
 covariatesfile       => 'data/covariates3std.txt',
 outcomevarfile       => 'data/outcomevars.txt',

 ##main options
 resultsdir       => 'results',
 displaylevel     => 3, #verbose output
 targetindicator  => 0, # diabetes in column 0
 #globalrho => 0,
 #outcomes => 1,
 samples  => 25,
 burnin   => 5,
 every    => 1,
 #indadmixhiermodel => 0,
 randommatingmodel  => 0,
 #fixedallelefreqs  => 1,
 numannealedruns    => 0,
 thermo => 0,

 ##output files
 logfile               => 'logfile.txt',
 paramfile             => 'paramfile.txt',
 regparamfile          => 'regparamfile.txt',
 indadmixturefile      => 'indadmixture.txt',
 ergodicaveragefile    => 'ergodicaverage.txt',
 #allelefreqoutputfile => 'allelefreqoutputfile.txt',

 ##optional tests
 #dispersiontestfile            => 'dispersiontest.txt',
 #admixturescorefile            => 'admixscorefile.txt',
 #residualallelicassocscorefile => 'resallelicassocscores.txt',
 allelicassociationscorefile    => 'allelicassociationscorefile.txt',
 #ancestryassociationscorefile  => 'ancestryassociationscorefile.txt',
 #affectedsonlyscorefile        => 'affectedsonlyscorefile.txt',
 haplotypeassociationscorefile  => 'hapassocscore.txt',
 stratificationtestfile         => 'strat_test.txt'
};

doAnalysis($prog, $rscript, $arg_hash);




