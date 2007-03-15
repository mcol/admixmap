#!/usr/bin/perl
# script to test most ADMIXMAP options and compare results with previous results
 
print "OS is ";print $^O;
$resultsdir = "results";

# Change this to the location of the admixmap executable
my $executable = './admixmap';

# $arg_hash is a hash of parameters passed to
# the executable as arguments.
#
# keys (left-hand side) are parameter names
# values (right-hand side) are parameter values
my $arg_hash = {
    samples                    => 15, 
    burnin                     => 5,
    every                      => 1,
    locusfile                  => 'data/loci.txt',
    genotypesfile              => 'data/genotypes.txt',
    outcomevarfile             => 'data/outcomevars.txt',
    covariatesfile             => 'data/covariates3std.txt',
    targetindicator            => 0, # diabetes in column 1
    displaylevel             => 2,
    
# output files
    logfile                    => 'logfile.txt',
    paramfile                  => 'param.txt',
    regparamfile               => 'regparam.txt',
    indadmixturefile           => 'indadmixture.txt',
    ergodicaveragefile         => 'ergodicaverage.txt',

# extra output files
    haplotypeassociationscorefile  => 'hapassocscore.txt',
    allelicassociationscorefile    => 'allelicassocscore.txt',
    allelefreqoutputfile           => 'allelefreqoutput.txt'
};

# single population, thermodynamic  
$arg_hash->{thermo} = 1;
$arg_hash->{numannealedruns} = 4;
$arg_hash->{populations} = 1;
#if(doAnalysis($executable,$arg_hash, $resultsdir)){
#    &CompareThenMove("results", "results0");
#}

# single population, reference prior on allele freqs, annealing  
$arg_hash->{thermo} = 0;
$arg_hash->{indadmixhiermodel} = 1;
$arg_hash->{stratificationtestfile}  = 'strat_test.txt';
#if(doAnalysis($executable,$arg_hash, $resultsdir)){
#  &CompareThenMove("results", "results1");
#}

# two populations, reference prior on allele freqs  
$arg_hash->{numannealedruns} = 0;
$arg_hash->{populations}     = 2;
#if(doAnalysis($executable,$arg_hash)){
#    &CompareThenMove("results", "results2");
#}

# fixed allele freqs, individual sumintensities
# possible problem here - changing numannealedruns should not change output except for annealmon.txt
delete $arg_hash->{populations};
delete $arg_hash->{allelefreqoutputfile};
$arg_hash->{fixedallelefreqs} = 1;
$arg_hash->{priorallelefreqfile}  = 'data/priorallelefreqs.txt',
$arg_hash->{allelefreqscorefile}  = 'allelefreqscorefile.txt';
$arg_hash->{allelefreqscorefile2} = 'allelefreqscorefile2.txt';
$arg_hash->{ancestryassociationscorefile} = 'ancestryassocscorefile.txt';
$arg_hash->{affectedsonlyscorefile}       = 'affectedsonlyscorefile.txt';
$arg_hash->{globalrho} = 0;
$arg_hash->{numannealedruns} = 5;
$arg_hash->{thermo} = 0;
#if(doAnalysis($executable,$arg_hash)){
    #&CompareThenMove("results", "results3");
#}

# prior on allele freqs
$arg_hash->{testoneindiv} = 0;
$arg_hash->{numannealedruns} = 0;
$arg_hash->{thermo} = 0;
$arg_hash->{fixedallelefreqs} = 0;
$arg_hash->{globalrho}        = 1;
$arg_hash->{randommatingmodel} = 1;
$arg_hash->{allelefreqoutputfile}     = 'allelefreqoutput.txt';
delete $arg_hash->{allelefreqscorefile};
delete $arg_hash->{allelefreqscorefile2};
delete $arg_hash->{affectedsonlyscorefile};
$arg_hash->{dispersiontestfile}  = 'dispersiontest.txt';
$arg_hash->{targetindicator} = 1; # skin reflectance
#if(doAnalysis($executable,$arg_hash)){
#    &CompareThenMove("results", "results4");
#}

# dispersion model for allele freqs 
delete $arg_hash->{priorallelefreqfile};
delete $arg_hash->{dispersiontestfile};
$arg_hash->{historicallelefreqfile} = 'data/priorallelefreqs.txt';
$arg_hash->{affectedsonlyscorefile}       = 'affectedsonlyscorefile.txt';
$arg_hash->{fstoutputfile} = 'FSToutputfile.txt';
$arg_hash->{dispparamfile} = 'disppar.txt';
$arg_hash->{randommatingmodel} = 0;
$arg_hash->{targetindicator} = 0; # diabetes
#if(doAnalysis($executable,$arg_hash)){
#    &CompareThenMove("results", "results5");
#}

# prior on allele freqs, testoneindiv, no regression, thermo
# possible problem here - changing numannealedruns should not change output except for annealmon.txt
delete $arg_hash->{historicallelefreqfile};
delete $arg_hash->{affectedsonlyscorefile};
delete $arg_hash->{ancestryassociationscorefile};
delete $arg_hash->{allelicassociationscorefile};
delete $arg_hash->{haplotypeassociationscorefile};
delete $arg_hash->{fstoutputfile};
delete $arg_hash->{dispparamfile};
delete $arg_hash->{regparamfile};
delete $arg_hash->{outcomes};
delete $arg_hash->{outcomevarfile};
delete $arg_hash->{covariatesfile};
$arg_hash->{testoneindiv} = 1;
$arg_hash->{priorallelefreqfile}  = 'data/priorallelefreqs.txt',
$arg_hash->{randommatingmodel} = 1;
$arg_hash->{globalrho} = 0;
$arg_hash->{testoneindiv} = 1;
$arg_hash->{numannealedruns} = 100;
$arg_hash->{thermo} = 1;
#if(doAnalysis($executable,$arg_hash)){
#    &CompareThenMove("results", "results6");
#}

# autosomal and X chromosome data
my $arg_hash = {
    genotypesfile                   => 'Cattledata/ngu_complete2.txt',
    locusfile                          => 'Cattledata/loci.txt',
    populations => 2,
    globalrho => 1,
    samples  => 25,
    burnin   => 5,
    every    => 1,
    resultsdir               => "$resultsdir",
    logfile                     => 'log.txt',
    paramfile               => 'paramfile.txt',
    indadmixturefile     => 'indadmixture.txt',
    ergodicaveragefile => 'ergodicaverages.txt',
    allelefreqoutputfile  => 'allelefreqs.txt',
    #allelicassociationscorefile       => 'allelicassociationscorefile.txt',
    #ancestryassociationscorefile  => 'ancestryassociationscorefile.txt',
    #affectedsonlyscorefile             => 'affectedsonlyscorefile.txt',
    #haplotypeassociationscorefile => 'hapassocscore.txt',
    #stratificationtestfile                   => 'strat_test.txt'
};
#if(doAnalysis($executable,$arg_hash)){
#    &CompareThenMove("results", "cattleresults");
#}

# Single individual
my $arg_hash = {
    burnin   => 10,
    samples  => 60,
    every    => 1,
    numannealedruns => 0, 
    displaylevel   => 2,

    locusfile                    => "IndData/loci.txt",
    genotypesfile                => "IndData/genotypes.txt",
    priorallelefreqfile          => "IndData/priorallelefreqs3way.txt",
    randommatingmodel            => 1,
    globalrho                    => 0,
    fixedallelefreqs             => 1,
    admixtureprior                   => "1,1,0",
    admixtureprior1                   => "1,1,1",
    logfile                      => "logfile.txt",
    chib                         => 1,
    indadmixturefile             => "indadmixture.txt"
};
#if(doAnalysis($executable,$arg_hash)){
#    &CompareThenMove("results", "Indresults");
#}

#hapmixmodel
my $executable = './hapmixmap2';
my $arg_hash = {
#data files

   genotypesfile                   => 'hapmixdata/genotypes_diploid.txt', #diploid data
#    genotypesfile                   => 'hapmixdata/genotypes_haploid.txt',#haploid data

    locusfile                       => 'hapmixdata/loci.txt',
    #priorallelefreqfile             => 'hapmixdata/allelefreqs.txt',
    #fixedallelefreqs                => 1,

    states=>4,
#    checkdata=> 0,

#main options
    resultsdir => 'results',
    displaylevel   => 3, 

    samples  => 15,
    burnin   => 5,
    every    => 1,

numannealedruns => 0,
thermo => 0,
hapmixmodel => 1,

hapmixlambdaprior=>"4000, 10, 100, 10",
#fixedmixtureprops=>1,

allelefreqprior => "1, 1, 1",
#initiallambdafile => "data/initialambdas.txt",
#allelefreqfile => "data/initialallelefreqs.txt",

lambdasamplerparams => "0.1, 0.00001, 10, 0.9, 20",

#output files
    logfile                     => 'logfile.txt',
    paramfile               => 'paramfile.txt',
    #regparamfile          => 'regparamfile.txt',
    #ergodicaveragefile => 'ergodicaverage.txt',

    finalallelefreqfile  => "finalallelefreqs.txt",
    finallambdafile => "finallambdas.txt",
    finalfreqpriorfile => "finalfreqpriors.txt",

    hapmixlambdaoutputfile =>"lambdaPosteriorMeans.txt",
    allelefreqprioroutputfile =>"etaPosteriorMeans.txt",

#optional tests
#residualallelicassocscorefile => 'residualLDscores.txt',
mhscoretestfile => 'mhtest.txt',
    #allelicassociationscorefile       => 'allelicassociationscorefile.txt',
};
#diploid data
if(doAnalysis($executable,$arg_hash)){
    &CompareThenMove("results", "HapmixDiploidresults");
}
#haploid data
$arg_hash->{genotypesfile} = "hapmixdata/genotypes_haploid.txt";
if(doAnalysis($executable,$arg_hash)){
    &CompareThenMove("results", "Hapmixresults");
}

#case-control analysis
$arg_hash->{ccgenotypesfile} = "hapmixdata/genotypes_casectrl.txt";
$arg_hash->{outcomevarfile} = "hapmixdata/outcome.txt";
$arg_hash->{allelicassociationscorefile} = "AllelicAssocScores.txt";
if(doAnalysis($executable,$arg_hash)){
    &CompareThenMove("results", "HapmixCCresults");
}


######################################################################################
sub doAnalysis {
    my ($prog, $args) = @_;
    my $command = $prog.getArguments($args);
    my$ status = system("$command");
    return status;
}

sub getArguments
{
    my $hash = $_[0];
#    my $arg = '';
#    foreach my $key (keys %$hash){
#	$arg .= ' --'. $key .'='. $hash->{$key};
#    }
#    return $arg;
    my $filename = 'perlargs.txt';
    open(OPTIONFILE, ">$filename") or die ("Could not open args file");
    foreach my $key (keys %$hash){
      print OPTIONFILE $key . '=' . $hash->{$key} . "\n";
    }
    close OPTIONFILE;
    return " ".$filename;
}

#SUBROUTINE TO COMPARE ALL FILES IN sourcedir WITH ORIGINALS AND MOVE
sub CompareThenMove {
    my ($sourcedir, $targetdir) = @_;
    my $prefix = "old_";
# define commands for different OS's
    if($^O eq "MSWin32") {$diffcmd = "fc"; $movecmd = "move /y"; $slash="\\";$delcmd="rmdir /s /q";}
    else {$diffcmd = "diff -s"; $movecmd = "mv -f"; $slash="/";$delcmd="rm -r";}
    if (-e $targetdir) { # compare with sourcedir
	opendir(SOURCE, $sourcedir) or die "can't open $sourcedir folder: $!";
	while ( defined (my $file = readdir SOURCE) ) {
	    next if (($file =~ /^\.\.?$/) || ($file eq "logfile.txt"));     # skip . and .. and logfile
	    system("$diffcmd $sourcedir$slash$file $targetdir$slash$prefix$file");
	    if($^O eq "MSWin32"){system("pause")}  # for checking comparisons 
	} #compare
	closedir(SOURCE);
    }
    if (-e $sourcedir) {
	if (-e $targetdir) {
	    my $olddir = "$prefix$targetdir";
	    if (-e $olddir){system("$delcmd $prefix$targetdir");} #delete old results directory (necessary for next command)
	    system("$movecmd $targetdir $prefix$targetdir"); #preserve old results by renaming directory
	    }  
	system("$movecmd $sourcedir $targetdir"); #rename results dir
	mkdir("$sourcedir");                      #we need results dir for next analysis
	opendir(TARGET, $targetdir) or die "can't open $targetdir folder: $!";
	while ( defined (my $file = readdir TARGET) ) { #for each results file
	    next if $file =~ /^\.\.?$/;     # skip . and ..
	    system("$movecmd $targetdir$slash$file $targetdir$slash$prefix$file ");} #prefix with "old_"
	closedir(TARGET);
    }
    else {print "$sourcedir does not exist\n"}
}
