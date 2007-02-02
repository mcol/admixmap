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

##default values for script: 4 states, 50 its with 10 burnin
my $samples=50;
my $burnin=10;
my $every=1;
my $STATES = 8;
my $POP = "Eur";
my $initial = 0;
my $cc = 0;

# Change this to the location of the admixmap executable
my $executable = './hapmixmap';

##parse any command line options
GetOptions(
	   "samples=i"=>\$samples, 
	   "burnin=i"=>\$burnin, 
	   "every=i"=>\$every, 
	   "states=i"=>\$STATES, 
	   "exec=s"=>\$executable, 
	   "init"=>\$initial, 
	   "cc"=>\$cc
);

#TODO: usage message and -h option

my $datadir = "data";
my $hapmapdir = "hapmap/$POP/data";

my $arg_hash = 
{
#main options
    resultsdir => 'results',
    displaylevel   => 3, 
    checkdata => 0,

    samples  => $samples,
    burnin   => $burnin,
    every    => $every,
    numannealedruns => 0,
    thermo   => 0,


#HapMap data files
    genotypesfile => "$hapmapdir/genotypes.txt",
    locusfile     => "$hapmapdir/loci.txt",

#model    
    #fixedallelefreqs => 1,
    states      => $STATES,
    hapmixmodel => 1,

#prior
    hapmixlambdaprior => "150, 1, 40, 10",    
    allelefreqprior   => "0.2, 1, 1",
    rhosamplerparams  => "0.1, 0.00001, 10, 0.9, 20",

#output files
    logfile            => 'logfile.txt',
    paramfile          => 'paramfile.txt',
    dispparamfile      => "allelefreqpriors.txt",
    #regparamfile      => 'regparamfile.txt',
    ergodicaveragefile => 'ergodicaverage.txt',

#saved state files
    allelefreqprioroutputfile => "initialetas.txt",
    allelefreqoutputfile      => "initialallelefreqs.txt",
    hapmixlambdaoutputfile    => "initiallambdas.txt",

#optional tests
};

##initial run
$arg_hash->{resultsdir}="Results$STATES"."States2";

#initial run
if($initial){
    doAnalysis($executable,$arg_hash);
}
else{
##rerun with final values of lambda, freqs in previous run as starting values
    if($cc){#running case-control data
	$arg_hash->{ccgenotypesfile} = "$datadir/genotypes.txt";
	$arg_hash->{outcomevarfile} = "$datadir/outcome.txt";
        $arg_hash->{allelicassociationscorefile} = 'allelicassocscores.txt';
    }
   $arg_hash->{initialfreqpriorfile} = "$datadir/initialetas.txt";
   $arg_hash->{initiallambdafile}  = "$datadir/initiallambdas.txt";
   $arg_hash->{initialallelefreqfile}           = "$datadir/initialallelefreqs.txt";
   doAnalysis($executable,$arg_hash);
}



print "script ended: ";
my $endtime = scalar(localtime());
print $endtime;
print "\n";



sub getArguments
{
    my $hash = $_[0];
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
   #my $command = "mpiexec " . $prog.getArguments($args);
    my $command = $prog.getArguments($args);

    $ENV{'RESULTSDIR'} = $args->{resultsdir};
    if(system($command)==0){
      system("cp $args->{resultsdir}/$args->{allelefreqoutputfile} $datadir");
      system("cp $args->{resultsdir}/$args->{hapmixlambdaoutputfile} $datadir");	
      system("cp $args->{resultsdir}/$args->{allelefreqprioroutputfile} $datadir");	
# Comment out the next three lines to run admixmap without R script
    print "Starting R script to process output\n";
      system("R CMD BATCH --quiet --no-save --no-restore ./AdmixmapOutput.R $args->{resultsdir}/Rlog.txt");
      print "R script completed\n\n";
  }
}



