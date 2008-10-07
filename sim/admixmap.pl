#!/usr/bin/perl -w
use strict;

my $function_file = "doanalysis.pl";

# Change this to the location of the admixmap executable
my $executable = 'admixmap';

# Change this to the location of the R script
my $rscript = "../dist/AdmixmapOutput.R";

##the following lines are a botch to make the script work straight out of the repository
if(!(-f $function_file) && (-f "../$function_file")){
##try one level up
    $function_file = "../$function_file";
}
sub doAnalysis
{
    my ($prog,$args) = @_;
    my $command = $prog.getArguments($args);
    if (-e $args->{resultsdir}) {
	rmtree($args->{resultsdir});
    } 
    mkpath($args->{resultsdir});
    $ENV{'RESULTSDIR'} = $args->{resultsdir};
    print "Results will be written to subdirectory $ENV{'RESULTSDIR'}";
    system($command);
    print "Starting R script to process output\n";
    system("R CMD BATCH --quiet --no-save --no-restore ../dist/AdmixmapOutput.R \
            $args->{resultsdir}/Rlog.txt");
      print "R script completed\n\n";
}
if(!(-f $rscript) && (-f "../$rscript")){
##try one level up
    $rscript = "../$rscript";
}
require $function_file or die("cannot find doanalysis.pl");

my $arg_hash = 
{
    samples                    => 250, 
    burnin                     => 50,
    every                      => 5,
    numannealedruns            => 0,
    locusfile                  => 'data/loci.txt',
    genotypesfile              => 'data/genotypes.txt',
    populations                => 2, 
    priorallelefreqfile        => 'data/priorallelefreqs.txt',
    outcomevarfile             => 'data/outcome.txt',
    displaylevel               => 3,
    logfile                    => 'log.txt',
    paramfile                  => 'param.txt',
    regparamfile               => 'regparam.txt',
    hwtest                     => 1,
};

# model with prior allele freqs
$arg_hash->{resultsdir}            = 'results';  
$arg_hash->{indadmixturefile}   = "indivadmixture.txt";
#$arg_hash->{dispersiontestfile}    = "dispersionTest.txt";
$arg_hash->{affectedsonlytest} = 1;
$arg_hash->{ancestryassociationtest} = 1;

print "script began: ";
my $starttime = scalar(localtime());
print $starttime;
print "\n";

doAnalysis($executable, $rscript, $arg_hash);

print "script ended: ";
my $endtime = scalar(localtime());
print $endtime;
print "\n";

