#!/usr/bin/perl
use strict; 
use File::Path;

sub getArguments
{
    my $hash = $_[0];
#    my $arg = '';
#   foreach my $key (keys %$hash){
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

sub doAnalysis
{
    my ($prog,$args) = @_;
    my $command = $prog.getArguments($args);
    
    $ENV{'RESULTSDIR'} = $args->{resultsdir};
    print "\nResults will be written to subdirectory $ENV{'RESULTSDIR'}\n";
    system($command);
    #print "Starting R script to process output\n";
    system("R CMD BATCH --quiet --no-save --no-restore ../test/AdmixmapOutput.R $args->{resultsdir}/Rlog.txt");
    #print "R script completed\n\n";
}

sub CompareThenMove {
    my ($sourcedir, $targetdir) = @_;
    my $prefix = "old_";
# define commands for different OS's
    my $diffcmd = "diff -s"; my $movecmd = "mv -f"; my $slash="/";my $delcmd="rm -R";
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

################### DO NOT EDIT ABOVE THIS LINE ########################
my $executable = '../test/admixmap';

my $arg_hash = {
#data files

#   genotypesfile                   => 'data/genotypes.txt', #diploid data
   genotypesfile                   => 'data/genotypes_haploid.txt',#haploid data

    locusfile                       => 'data/loci.txt',
    #priorallelefreqfile             => 'data/allelefreqs.txt',
    #fixedallelefreqs                => 1,

    populations=>6,

#main options
    resultsdir => 'results',
    displaylevel   => 3, 

    samples  => 150,
    burnin   => 50,
    every    => 5,

    numannealedruns => 0,
    thermo => 0,
    hapmixmodel => 1,
#   indadmixhiermodel => 0,
    randommatingmodel => 0,
    checkdata=> 0,

hapmixlambdaprior=>"400, 1, 10, 1",

allelefreqprior => "2, 10, 1",
#initialhapmixlambdafile => "data/initialambdas.txt",
#allelefreqfile => "data/initialallelefreqs.txt",

rhosamplerparams => "0.5, 0.00001, 10, 0.9, 20",

#output files
    logfile                     => 'logfile.txt',
    paramfile               => 'paramfile.txt',
    #regparamfile          => 'regparamfile.txt',
    #indadmixturefile     => 'indadmixture.txt',
    #ergodicaveragefile => 'ergodicaverage.txt',
    allelefreqoutputfile  => "initialallelefreqs.txt",
allelefreqoutputfile =>"allelefreqpriors.txt",
    hapmixlambdaoutputfile => "data/initiallambdas.txt",

#optional tests
residualallelicassocscorefile => 'residualLDscores.txt',
    #allelicassociationscorefile       => 'allelicassociationscorefile.txt',
};

# Initial run 
#haploid data
$arg_hash->{resultsdir}            = 'ResultsHaploid';  
doAnalysis($executable,$arg_hash);
system("cp Results/initialallelefreqs.txt data");
#CompareThenMove("Results", "Results4");

#diploid data
$arg_hash->{resultsdir}            = 'ResultsDiploid';  
$arg_hash->{genotypesfile} = "data/genotypes.txt";
doAnalysis($executable,$arg_hash);
#system("cp Results/initialallelefreqs.txt data");
#CompareThenMove("Results", "Results4");

# rerun with final values of previous run as intial values of this
$arg_hash->{allelefreqfile}="data/initialallelefreqs.txt";
$arg_hash->{initialhapmixlambdafile}="data/initiallambdas.txt";
$arg_hash->{fixedallelefreqs} = 0;
delete $arg_hash->{priorallelefreqfile};
#doAnalysis($executable,$arg_hash);
