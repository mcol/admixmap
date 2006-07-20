#!/usr/bin/perl -w
use strict;

# Change this to the location of the admixmap executable
my $executable = './admixmap';

# $arg_hash is a hash of parameters passed to
# the executable as arguments.
#
# keys (left-hand side) are parameter names
# values (right-hand side) are parameter values
my $arg_hash = 
{
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

print "script began: ";
my $starttime = scalar(localtime());
print $starttime;
print "\n";

doAnalysis($executable,$arg_hash);

print "script ended: ";
my $endtime = scalar(localtime());
print $endtime;
print "\n";

sub getArguments
{
    my $hash = $_[0];
    my $arg = '';
    foreach my $key (keys %$hash){
	$arg .= ' --'. $key .'='. $hash->{$key};
    }
    return $arg;
}

sub doAnalysis {
    my ($prog,$args) = @_;
    my $command = $prog.getArguments($args);

    $ENV{'RESULTSDIR'} = $args->{resultsdir};
    print "\nResults will be written to subdirectory $ENV{'RESULTSDIR'}\n";
    my $status = system($command);

    # Comment out the remaining lines to run admixmap without R script
    if( $status == 0)
      {
	my $rcmd = "R CMD";
	if($^O eq "MSWin32") {
	  $rcmd = "Rcmd";
	}
	print "Starting R script to process output\n";
	system("$rcmd BATCH --quiet --no-save --no-restore ../test/AdmixmapOutput.R $args->{resultsdir}/Rlog.txt\n");
	print "R script completed\n\n";
      }else{
	print "Warning: admixmap has not run successfully, R script will not be run.\n\n"
      }
}



