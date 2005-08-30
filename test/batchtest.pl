#!/usr/bin/perl
# Test Script to test most of ADMIXMAP options and compare results with previous results
# 
print "OS is ";print $^O;
$resultsdir = "results";
if (-e $resultsdir){ 
system("erase /q $resultsdir");
system("mkdir $resultsdir");}
else {system("mkdir $resultsdir");}

# Change this to the location of the admixmap executable
my $executable = 'admixmap';

# $arg_hash is a hash of parameters passed to
# the executable as arguments.
#
# keys (left-hand side) are parameter names
# values (right-hand side) are parameter values
my $arg_hash = 
{
    samples                    => 20, 
    burnin                     => 5,
    every                      => 1,
    locusfile                  => 'data/loci.txt',
    genotypesfile              => 'data/genotypes.txt',
    outcomevarfile             => 'data/outcomevars.txt',
    covariatesfile             => 'data/covariates3std.txt',
    targetindicator            => 0, # diabetes in column 1
    analysistypeindicator      => 5, # one binary and one continuous outcome var 
    coutindicator              => 1,
    
# output files
    logfile                    => 'logfile.txt',
    paramfile                  => 'param.txt',
    regparamfile  => 'regparam.txt',
    indadmixturefile           => 'indadmixture.txt',
    ergodicaveragefile         => 'ergodicaverage.txt',

# extra output files
    haplotypeassociationscorefile  => 'hapassocscore.txt',
    allelicassociationscorefile    => 'allelicassocscore.txt',
    stratificationtestfile         => 'strat_test.txt',
    allelefreqoutputfile           => 'allelefreqoutput.txt'
};

# single population, reference prior on allele freqs  
$arg_hash->{populations} = 1;
doAnalysis($executable,$arg_hash, $resultsdir);
&CompareThenMove("results", "results1");

# two populations, reference prior on allele freqs  
$arg_hash->{populations}     = 2;
$arg_hash->{sumintensitiesalpha}  = 5;
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
$arg_hash->{dispersiontestfile}  = 'dispersiontest.txt';
$arg_hash->{analysistypeindicator} = 2; # continuous outcome var
$arg_hash->{targetindicator} = 1; # skin reflectance
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
$arg_hash->{analysistypeindicator} = 3; # binary outcome var
$arg_hash->{targetindicator} = 0; # diabetes
doAnalysis($executable,$arg_hash);
&CompareThenMove("results", "results5");

#Single individual
my $arg_hash = 
{
    burnin   => 10,
    samples  => 51,
    every    => 2,
    analysistypeindicator     => -1,  
    targetindicator => 1, # offset (from column 1) of column containing outcome variable
    coutindicator   => 1,

    locusfile                    => "IndData/loci.txt",
    genotypesfile                => "IndData/genotypes.txt",
    priorallelefreqfile          => "IndData/priorallelefreqs3way.txt",
    randommatingmodel            => 1,
    globalrho                    => 0,
    sumintensitiesalpha => 1.0, #flat prior on sumintensities
    sumintensitiesbeta => 0.0,

#    fixedallelefreqs             => 1,
    initalpha0                   => "1,1,1",
    initalpha1                   => "1,1,0",

    logfile                      => "logfile.txt",
    marglikelihood => 1,
    indadmixturefile             => "indadmixture.txt"
};

 doAnalysis($executable,$arg_hash);
 &CompareThenMove("results", "Indresults");

sub doAnalysis
{
    my ($prog, $args) = @_;
    my $command = $prog.getArguments($args);
    system("$command");
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
#define commands for different OS's
    if($^O eq "MSWin32") {$diffcmd = "fc"; $movecmd = "move /y"; $slash="\\";$delcmd="rmdir /s /q";}
    else {$diffcmd = "diff -s"; $movecmd = "mv -f"; $slash="/";$delcmd="rm -r";}
##
    if (-e $targetdir) { # compare with sourcedir
	opendir(SOURCE, $sourcedir) or die "can't open $sourcedir folder: $!";
	while ( defined (my $file = readdir SOURCE) ) {
	    next if (($file =~ /^\.\.?$/) || ($file eq "logfile.txt"));     # skip . and .. and logfile
	    system("$diffcmd $sourcedir$slash$file $targetdir$slash$prefix$file");
	if($^O eq "MSWin32"){system("pause")}  # for checking comparisons 
	}#compare

	closedir(SOURCE);

	}
###################
    if (-e $sourcedir) {
         if (-e $targetdir) {
	 my $olddir = "$prefix$targetdir";
	 if (-e $olddir){system("$delcmd $prefix$targetdir");} #delete old results directory (necessary for next command)
         system("$movecmd $targetdir $prefix$targetdir");#preserve old results by renaming directory
	     }  

	 system("$movecmd $sourcedir $targetdir"); #rename results dir
	 mkdir("$sourcedir");                      #we need results dir for next analysis
	 opendir(TARGET, $targetdir) or die "can't open $targetdir folder: $!";
	 while ( defined (my $file = readdir TARGET) ) {#for each results file
	     next if $file =~ /^\.\.?$/;     # skip . and ..
	     system("$movecmd $targetdir$slash$file $targetdir$slash$prefix$file ");} #prefix with "old_"
	 closedir(TARGET);
     }
    else {print "$sourcedir does not exist\n"}
}
