#!/usr/bin/perl
# Script to test most ADMIXMAP options and compare with previous results
use strict;
use File::Path;
use File::Copy;
print "OS is ";print $^O;

# Change this to the location of the admixmap executable
my $executable = './admixmap';

my $arg_hash = {
    samples                    => 20, 
    burnin                     => 5,
    every                      => 1,
    locusfile                  => 'data/loci.txt',
    genotypesfile              => 'data/genotypes.txt',
    outcomevarfile             => 'data/outcomevars.txt',
    covariatesfile             => 'data/covariates3std.txt',
    targetindicator            => 0, # diabetes in column 1
    coutindicator              => 1,
    
# output files
    resultsdir                 => 'results',
    logfile                    => 'logfile.txt',
    paramfile                  => 'param.txt',
    regparamfile               => 'regparam.txt',
    indadmixturefile           => 'indadmixture.txt',
    ergodicaveragefile         => 'ergodicaverage.txt',

# extra output files
    haplotypeassociationscorefile  => 'hapassocscore.txt',
    allelicassociationscorefile    => 'allelicassocscore.txt',
    stratificationtestfile         => 'strat_test.txt'
};

# single population, reference prior on allele freqs  
$arg_hash->{populations} = 1;
doAnalysis($executable,$arg_hash);
CompareThenMove("results", "results1");

# two populations, reference prior on allele freqs  
$arg_hash->{populations}     = 2;
$arg_hash->{globalrho}      = 0; # should not be necessary if option sumintensitiesprior is specified with 3 elements
$arg_hash->{sumintensitiesprior}  = "8,11,10";
doAnalysis($executable,$arg_hash);
&CompareThenMove("results", "results2");

# fixed allele freqs
delete $arg_hash->{populations};
$arg_hash->{fixedallelefreqs} = 1;
$arg_hash->{priorallelefreqfile}  = 'data/priorallelefreqs.txt',
$arg_hash->{allelefreqscorefile}  = 'allelefreqscorefile.txt';
$arg_hash->{allelefreqscorefile2} = 'allelefreqscorefile2.txt';
$arg_hash->{ancestryassociationscorefile} = 'ancestryassocscorefile.txt';
$arg_hash->{affectedsonlyscorefile}       = 'affectedsonlyscorefile.txt';
$arg_hash->{globalrho} = 0;
doAnalysis($executable,$arg_hash);
&CompareThenMove("results", "results3");

# prior on allele freqs
$arg_hash->{fixedallelefreqs} = 0;
$arg_hash->{globalrho}        = 1;
$arg_hash->{randommatingmodel} = 1;
delete $arg_hash->{allelefreqscorefile};
delete $arg_hash->{allelefreqscorefile2};
delete $arg_hash->{affectedsonlyscorefile};
$arg_hash->{allelefreqoutputfile}  = 'allelefreqoutput.txt';
$arg_hash->{dispersiontestfile}  = 'dispersiontest.txt';
$arg_hash->{targetindicator} = 1; # skin reflectance
$arg_hash->{outcomes}=1;
doAnalysis($executable,$arg_hash);
&CompareThenMove("results", "results4");

# dispersion model for allele freqs 
delete $arg_hash->{priorallelefreqfile};
delete $arg_hash->{dispersiontestfile};
$arg_hash->{historicallelefreqfile} = 'data/priorallelefreqs.txt';
$arg_hash->{affectedsonlyscorefile}       = 'affectedsonlyscorefile.txt';
$arg_hash->{fstoutputfile} = 'FSToutputfile.txt';
$arg_hash->{dispparamfile} = 'disppar.txt';
$arg_hash->{randommatingmodel} = 0;
$arg_hash->{targetindicator} = 0; # diabetes
$arg_hash->{outcomes}=1;
doAnalysis($executable,$arg_hash);
&CompareThenMove("results", "results5");

#Single individual
my $arg_hash = {
    burnin   => 10,
    samples  => 51,
    every    => 2,
    coutindicator   => 1,
    locusfile                    => "IndData/loci.txt",
    genotypesfile                => "IndData/genotypes.txt",
    priorallelefreqfile          => "IndData/priorallelefreqs3way.txt",

    indadmixhiermodel => 0,
    randommatingmodel            => 1,
    globalrho                    => 0,
    sumintensitiesprior => "6,5,4", 

    initalpha0                   => "1,1,1",
    initalpha1                   => "1,1,0",

    resultsdir                   => 'results',
    logfile                      => "logfile.txt",
    marglikelihood => 1,
    indadmixturefile             => "indadmixture.txt"
};

 doAnalysis($executable,$arg_hash);
 &CompareThenMove("results", "Indresults");

sub doAnalysis
{
    my ($prog,$args) = @_;
    my $command = $prog.getArguments($args);
    print $args->{resultsdir};
    unless (-e "$args->{resultsdir}"){
	system("mkdir $args->{resultsdir}");
    }

    $ENV{'RESULTSDIR'} = $args->{resultsdir};
    print "Results will be written to subdirectory $ENV{'RESULTSDIR'}";
    system($command);
}

sub getArguments
{
    my $hash = $_[0];
    my $arg = '';
    foreach my $key (keys %$hash){
	$arg .= ' --'. $key .'='. $hash->{$key};
    }
    return $arg;
}

#SUBROUTINE TO COMPARE ALL FILES IN sourcedir WITH ORIGINALS 
# AND MOVE
sub CompareThenMove{
    my ($sourcedir, $targetdir) = @_;
    my $prefix = "old_";
    # define commands for different OS's
    # should redo this using perl functions that work with any OS 
    my $diffcmd;
    my $movecmd;
    my $slash;
    my $delcmd;
    if($^O eq "MSWin32") {$diffcmd = "fc"; $movecmd = "move /y"; $slash="\\"; $delcmd="rmdir /s /q";}
    else {$diffcmd = "diff -s"; $movecmd = "mv -f"; $slash="/"; $delcmd="rm -r";}
    if (-e $targetdir) { # compare with sourcedir
	opendir(SOURCE, $sourcedir) or die "can't open $sourcedir folder: $!";
	while ( defined (my $file = readdir SOURCE) ) {
	    next if (($file =~ /^\.\.?$/) || ($file eq "logfile.txt"));     # skip . and .. and logfile
	    system("$diffcmd $sourcedir$slash$file $targetdir$slash$prefix$file");
	    if($^O eq "MSWin32"){system("pause")}  # for checking comparisons 
	} # compare
	    closedir(SOURCE);
    }
    if (-e $sourcedir) { # if results directory exists
	if (-e $targetdir) { # if resultsn directory exists
	    my $olddir = "$prefix$targetdir"; # = old_resultsn
	    if (-e $olddir) { # remove old_resultsn directory
		rmtree $olddir;
	    } 
	    rename $targetdir, $olddir or die("Move from $targetdir to $olddir failed"); 
	} 
	rename $sourcedir, $targetdir or die("Move from $sourcedir to $targetdir failed");
	mkdir("$sourcedir");                      
	opendir(TARGET, $targetdir) or die("Can't open $targetdir folder: $!");
	while ( defined (my $file = readdir TARGET) ) {
	    ## for each results file
	    next if $file =~ /^\.\.?$/;     # skip . and ..
	    system("$movecmd $targetdir$slash$file $targetdir$slash$prefix$file ");} #prefix with "old_"
	closedir(TARGET);
    }
    else {print "$sourcedir does not exist\n"}
}
