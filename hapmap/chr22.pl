#!/usr/bin/perl -w
use strict;
use File::Path;
use Time::Local;

my $DEBUG = 0; # zero gives less output
print "script began: ";
my $starttime = scalar(localtime());
print $starttime;
print "\n";

if($#ARGV < 1){
    print("Invalid args\n");
    exit;
}
my $POP;
if($ARGV[0]==1){$POP="Eur";}
else {
    if($ARGV[0]==2){$POP="Afr";}
    else {if($ARGV[0]==3){$POP="Asian";}
	  else{
	      print ("Invalid population code - must be 1, 2 or 3\n");
	      exit;
	  }
      }
}
my $STATES = $ARGV[1];

# Change this to the location of the admixmap executable
#my $executable = '../test/adm-para';
#my $executable = '../test/admixmap';
my $executable = '/ichec/home/users/doducd/test/admixmap';


# $arg_hash is a hash of parameters passed to
# the executable as arguments.
#
# keys (left-hand side) are parameter names
# values (right-hand side) are parameter values
my $datadir = "/ichec/work/ndlif006b/genepi/hapmap/$POP/chr22data";
my $arg_hash = 
{
#data files
    genotypesfile                   => "$datadir/genotypes.txt",
    locusfile                          => "$datadir/loci.txt",
    #priorallelefreqfile             => 'data/priorallelefreqs.txt',
    #fixedallelefreqs => 1,
    populations=>$STATES,
    #outcomevarfile => 'chr22/dummyoutcome.txt',

checkdata=> 0,

#main options
    resultsdir => 'results',
    displaylevel   => 3, 

    samples  => 250,
    burnin   => 50,
    every    => 1,

numannealedruns => 0,
thermo => 0,
hapmixmodel => 1,
#indadmixhiermodel => 0,
randommatingmodel => 0,

hapmixlambdaprior=>"4, 0.1, 1, 1",

allelefreqprior => "1, 1, 1",
#initialhapmixlambdafile => "$datadir/initialambdas.txt",
#allelefreqfile => "$datadir/initialallelefreqs.txt",

rhosamplerparams => "0.1, 0.00001, 10, 0.9, 20",

#output files
    logfile                     => 'logfile.txt',
    paramfile               => 'paramfile.txt',
    #regparamfile          => 'regparamfile.txt',
    #indadmixturefile     => 'indadmixture.txt',
    #ergodicaveragefile => 'ergodicaverage.txt',
    allelefreqoutputfile  => "initialallelefreqs.txt",
    hapmixlambdaoutputfile => "$datadir/initiallambdas.txt",

#optional tests
#residualallelicassocscorefile => 'residualLDscores.txt',
    #allelicassociationscorefile       => 'allelicassociationscorefile.txt',
};

#model with $STATES block states

##initial run
$arg_hash->{resultsdir}="$POP/Results$STATES"."States";
doAnalysis($executable,$arg_hash);

##rerun with final values of lambda, freqs in previous run as starting values
$arg_hash->{resultsdir}="$POP/Results$STATES"."States2";
$arg_hash->{initialhapmixlambdafile} = "$datadir/initiallambdas.txt";
$arg_hash->{allelefreqfile} = "$datadir/initialallelefreqs.txt";
$arg_hash->{residualallelicassocscorefile} = 'residualLDscores.txt';
doAnalysis($executable,$arg_hash);

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


## half genotypes missing at half the loci
$arg_hash->{resultsdir}='MissingResults';
$arg_hash->{genotypesfile} = 'MissingData/genotypes.txt';
$arg_hash->{locusfile}  = 'MissingData/loci.txt';
$arg_hash->{outcomevarfile} = 'smallMissingData/outcome.txt';# dummy outcome for allelic assoc test
#$arg_hash->{allelicassociationscorefile} = 'allelicassociationscorefile.txt';
#doAnalysis($executable,$arg_hash);

## first 1000 loci, half genotypes missing at half the loci
$arg_hash->{resultsdir}='smallMissingResults';
$arg_hash->{genotypesfile} = 'smallMissingData/genotypes.txt';
$arg_hash->{locusfile}  = 'smallMissingData/loci.txt';
#doAnalysis($executable,$arg_hash);

print "script ended: ";
my $endtime = scalar(localtime());
print $endtime;
print "\n";



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

sub doAnalysis
{
    my ($prog,$args) = @_;
   #my $command = "mpiexec " . $prog.getArguments($args);
    my $command = $prog.getArguments($args);

    $ENV{'RESULTSDIR'} = $args->{resultsdir};
    if(system($command)==0){
	if($arg_hash->{allelefreqoutputfile}){
#copy file with final lambdas to data dir ready for next time
	system("cp $arg_hash->{resultsdir}/$arg_hash->{allelefreqoutputfile} $datadir/$arg_hash->{allelefreqoutputfile}");
    }
# Comment out the next three lines to run admixmap without R script
    print "Starting R script to process output\n";
    system("R --quiet --no-save --no-restore <../test/AdmixmapOutput.R >$args->{resultsdir}/Rlog.txt RESULTSDIR=$args->{resultsdir}");
    print "R script completed\n\n";
}
}



