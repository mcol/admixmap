#!/usr/bin/perl -w
use strict;

my $executable = './admixmap';
my $arg_hash = {
#data files
    genotypesfile                   => 'simdata/genotypes.txt',
    locusfile                       => 'simdata/loci.txt',
    #priorallelefreqfile             => 'simdata/priorallelefreqs.txt',
    # covariatesfile                  => 'data/covariates3.txt',
    outcomevarfile                  => 'simdata/outcome.txt',
fixedallelefreqs =>1,
#main options
    displaylevel   => 3, #verbose output
    samples  => 250,
    burnin   => 50,
    every    => 5,
    populations => 2,
    numannealedruns => 0,
    #indadmixhiermodel => 0,
    #admixtureprior => "3,1",

#output files
    resultsdir               => "simResults",
    logfile                  => 'logfile.txt',
    paramfile                => 'paramfile.txt',
    #regparamfile             => 'regparam.txt',
    indadmixturefile         => 'indadmixture.txt',
    #ergodicaveragefile       => 'ergodicaverage.txt',
    allelefreqoutputfile  => 'allelefreqoutputfile.txt',

#optional tests
     allelicassociationscorefile       => 'allelicassociationscorefile.txt',
     #ancestryassociationscorefile  => 'ancestryassociationscorefile.txt',
    #affectedsonlyscorefile             => 'affectedsonlyscorefile.txt',
     #haplotypeassociationscorefile => 'hapassocscore.txt',
    dispersiontestfile                   => 'dispersiontest.txt',
    #hwscoretestfile                  => 'HardyWeinbergTest.txt',
    stratificationtestfile           => 'stratificationtest.txt',
    #residualallelicassocscorefile    => 'residualLDscoretest.txt'
};

doAnalysis($executable,$arg_hash);

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

sub doAnalysis {
    my ($prog,$args) = @_;
    my $command = $prog.getArguments($args);
    $ENV{'RESULTSDIR'} = $args->{resultsdir};
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
