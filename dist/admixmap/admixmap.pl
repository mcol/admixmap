#!/usr/bin/perl -w
use strict;

my $function_file = "doanalysis.pl";

# Change this to the location of the admixmap executable
my $executable = 'admixmap';

# Change this to the location of the R script
my $rscript = "AdmixmapOutput.R";

##the following lines are a botch to make the script work straight out of the repository
if(!(-f $function_file) && (-f "../$function_file")){
##try one level up
  $function_file = "../$function_file";
}
if(!(-f $rscript) && (-f "../$rscript")){
##try one level up
  $rscript = "../$rscript";
}

require $function_file or die("cannot find doanalysis.pl");

# $arg_hash is a hash of parameters passed to
# the executable as arguments.
#
# keys (left-hand side) are parameter names
# values (right-hand side) are parameter values
my $arg_hash = 
{
 ##data files
 genotypesfile        => 'tutorial/data/genotypes.txt',
 locusfile            => 'tutorial/data/loci.txt',
 priorallelefreqfile  => 'tutorial/data/priorallelefreqs.txt',
 #populations => 1
 covariatesfile       => 'tutorial/data/covariates3std.txt',
 outcomevarfile       => 'tutorial/data/outcomevars.txt',

 ##main options
 resultsdir       => 'results',
 displaylevel     => 3, #verbose output
 #globalrho => 0,
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
 #dispersiontest            => 1,
 #admixtureassoctest        => 1,
 #residualldtest            => 1,
 allelicassociationtest     => 1,
 #ancestryassociationtest   => 1,
 #affectedsonlytest         => 1,
 haplotypeassociationtest   => 1,
 stratificationtest         => 1
};

print "script began: ";
my $starttime = scalar(localtime());
print $starttime;
print "\n";

doAnalysis($executable, $rscript, $arg_hash);

print "script ended: ";
my $endtime = scalar(localtime());
print $endtime;
print "\n";




