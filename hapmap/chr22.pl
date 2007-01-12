#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Path;
use Time::Local;

my $DEBUG = 0; # zero gives less output
print "script began: ";
my $starttime = scalar(localtime());
print $starttime;
print "\n";

##default values for script: Europeans, 4 states, 50 its with 10 burnin
my $samples=50;
my $burnin=10;
my $every=1;
my $POP = "Eur";
my $STATES = 4;
my $parallel = '';
my $maskfile = '';

# Change this to the location of the admixmap executable
#my $executable = '../test/adm-para';

#my $executable = '../test/hapmixmap';
my $executable = '/ichec/home/users/doducd/test/hapmixmap';

##parse any command line options
GetOptions("parallel" =>\$parallel, "samples=i"=>\$samples, "burnin=i"=>\$burnin, "every=i"=>\$every, "pop=s"=>\$POP, "states=i"=>\$STATES, "exec=s"=>\$executable, "maskfile=s" =>\$maskfile);

my $datadir = "$POP/chr22data";
# $arg_hash is a hash of parameters passed to
# the executable as arguments.
#
# keys (left-hand side) are parameter names
# values (right-hand side) are parameter values
my $arg_hash = 
{
#data files
#    genotypesfile                   => "$datadir/genotypes5000.txt",
#    locusfile                          => "$datadir/loci5000.txt",

#phased data
    genotypesfile                   => "$datadir/genotypes_phased.txt",
    locusfile                          => "$datadir/loci_phased.txt",

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

#model with $STATES block states

##initial run
$arg_hash->{resultsdir}="$POP/Results2$STATES";
doAnalysis($executable,$arg_hash);

##rerun with final values of lambda, freqs in previous run as starting values
$arg_hash->{resultsdir}="$POP/Results$STATES"."States";
$arg_hash->{initialhapmixlambdafile} = "$datadir/initiallambdas.txt";
$arg_hash->{allelefreqfile} = "$datadir/initialallelefreqs.txt";
$arg_hash->{residualallelicassocscorefile} = 'residualLDscores.txt';
#doAnalysis($executable,$arg_hash);

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

print "script ended: ";
my $endtime = scalar(localtime());
print $endtime;
print "\n";



sub getArguments
{
    my $hash = $_[0];
    my $filename = "args$POP$STATES.txt";
    open(OPTIONFILE, ">$filename") or die ("Could not open args file");
    foreach my $key (keys %$hash){
      print OPTIONFILE $key . '=' . $hash->{$key} . "\n";
    }

    if($maskfile){##if we are using masked data
	# If possible, append the contents of the index file with masked
	# loci to the option file, so users don't have to do it by hand.
	my $line;
	open(EXTERNAL_ARGS, "$datadir/$maskfile")
	    or warn("Can't open the external arguments file.");
	foreach $line (<EXTERNAL_ARGS>) {
	    print OPTIONFILE $line;
	}
	close(EXTERNAL_ARGS);
    }
    close OPTIONFILE;
    return " ".$filename;
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
	
# Comment out the next three lines to run admixmap without R script
    print "Starting R script to process output\n";
    system("R CMD BATCH --quiet --no-save --no-restore ../test/AdmixmapOutput.R $args->{resultsdir}/Rlog.txt RESULTSDIR=$args->{resultsdir}");
    print "R script completed\n\n";
}
}
