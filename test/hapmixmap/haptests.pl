#!/usr/bin/perl
use strict; 
use File::Path;
use Getopt::Long;

my $parallel = '';
my $simulate = '';
my $executable = '';

my $serial_executable = "$ENV{'HOME'}/usr/bin/hapmixmap";
my $parallel_executable = "$ENV{'HOME'}/usr/bin/hapmixmap-para";
##NB make sure this matches the name in the R script
my $datadir = "batchtestdata";

my $usage = 0;

GetOptions("exec=s"    => \$executable,
           "help!"     => \$usage,
           "parallel"  => \$parallel,
	   "simulate"  => \$simulate);

if ($usage) {
    print "\n";
    print "Usage: $0 [arguments]\n";
    print "  --exec <file>              Name of the HAPMIXMAP executable,\n";
    print "                             e.g. ../test/hapmixmap\n";
    print "  --parallel                 run parallel version of HAPMIXMAP,\n";
    print "                             when on a cluster or multi-processor.\n";
    print "  --simulate                 Simulate data first\n";
    print "\n";
    exit(1);
}

if($simulate){
##write settings file for R script
  open(RINPUTFILE, ">simHapMixSettings.R") or die("Could not open Rinputfile file");
  print RINPUTFILE "datadir<-\"$datadir\"\n";
  close(RINPUTFILE);

  print "Running R script to simulate data\n";
  if(system("R CMD BATCH --vanilla simHapMix.R Rsimlog.txt")){
    die "R script failed";
  }else{
    print "simulation complete\n";
  }
}

if(!$executable){
  if($parallel){
    $executable = $parallel_executable;
  }else{
    $executable = $serial_executable;
  }
}
my $function_file = "../../dist/doanalysis.pl";

require $function_file or die("cannot find doanalysis.pl");

my $copy = "cp";
my $slash="/";
if($^O eq "MSWin32"){
  $copy = "copy";
  $slash="\\";
}

my $arg_hash = 
{
 #data files
 locusfile                       => "$datadir/loci.txt",
 genotypesfile                   => "$datadir/genotypes_haploid.txt",#haploid data

 #main options
 displaylevel   => 3,
 samples  => 25,
 burnin   => 5,
 every    => 1,
 numannealedruns => 0,
 thermo => 0,
 checkdata=> 0,
 
 #model
 hapmixmodel => 1,
 states=>8,
 
 #priors
 arrivalrateprior=>"400, 1, 10, 1",
 allelefreqprecisionprior => "2, 10, 1",
 freqprecisionhiermodel => 1,
 
 arrivalratesamplerparams => "0.5, 0.00001, 10, 0.9, 20",
 
 #output files
 resultsdir => 'results',
 logfile                     => 'logfile.txt',
 paramfile               => 'paramfile.txt',
 freqprecisionfile => "allelefreqpriorsamples.txt",
 #regparamfile          => 'regparamfile.txt',
 ergodicaveragefile => 'ergodicaverage.txt',
 
#final values
 finalallelefreqfile  => "initialallelefreqs.txt",
 finalfreqpriorfile   =>"initialallelefreqpriors.txt",
 finalarrivalratefile =>"initialarrivalratess.txt",
 finalmixturepropsfile => "initialmixtureprops.txt",
 
 #posterior means
 arrivalrateposteriormeanfile => "lambdaPosteriorMeans.txt",
 allelefreqprecisionposteriormeanfile => "freqDispersionPosteriorMeans.txt",
 
 #optional tests
 mhscoretestfile => 'MantelHaenszelTest.txt'
 #residualallelicassocscorefile => 'residualLDscores.txt',
 #allelicassociationscorefile       => 'allelicassociationscorefile.txt',
};

#haploid data, from scratch, fixed allele freqs
$arg_hash->{resultsdir}            = 'ResultsHaploidFixedFreqs';
$arg_hash->{priorallelefreqfile}   = "$datadir/allelefreqs.txt";
$arg_hash->{fixedallelefreqs} = 1;
$arg_hash->{fixedmixturepropsprecision} = 0;
callDoAnalysis();
CompareThenMove($arg_hash->{resultsdir});


#diploid data, random allele freqs with default prior, fixed mixture props
$arg_hash->{resultsdir}            = 'ResultsDiploidFixedMixtureProps';
$arg_hash->{genotypesfile} = "$datadir/genotypes_diploid.txt";
$arg_hash->{fixedmixtureprops} = 1;
$arg_hash->{fixedmixturepropsprecision} = 1;
$arg_hash->{fixedallelefreqs} = 0;
delete $arg_hash->{priorallelefreqfile};
callDoAnalysis();
CompareThenMove($arg_hash->{resultsdir});

#haploid, random allele freqs with default prior, random mixture props
$arg_hash->{genotypesfile} = "$datadir/genotypes_haploid.txt";
$arg_hash->{resultsdir}            = 'ResultsHaploid';
$arg_hash->{fixedmixtureprops} = 0;
callDoAnalysis();
CompareThenMove($arg_hash->{resultsdir});

#haploid, resume
#system("$copy $arg_hash->{resultsdir}$slash"."initial*.txt $datadir");

$arg_hash->{resultsdir}            = 'ResultsHaploidRerun';
$arg_hash->{initialallelefreqfile}="$datadir/initialallelefreqs.txt";
$arg_hash->{initialarrivalratefile}="$datadir/initialarrivalrates.txt";
$arg_hash->{initialfreqpriorfile} = "$datadir/initialfreqpriors.txt";
$arg_hash->{initialarrivalratefile}="$datadir/initialmixtureprops.txt";
callDoAnalysis();
CompareThenMove($arg_hash->{resultsdir});


#Case-control analysis
delete $arg_hash->{initialallelefreqfile};
delete $arg_hash->{initialarrivalratefile};
delete $arg_hash->{initialfreqpriorfile};
delete $arg_hash->{initialarrivalratefile};
$arg_hash->{resultsdir}      = 'ResultsCaseControl';
$arg_hash->{ccgenotypesfile} = "$datadir/genotypes_casectrl.txt";
$arg_hash->{outcomevarfile} = "$datadir/outcome.txt";
$arg_hash->{allelicassociationscorefile} = 'AllelicAssocTests.txt';
callDoAnalysis();
CompareThenMove($arg_hash->{resultsdir});

sub callDoAnalysis {
    if($parallel){
	##doParallelAnalysis($executable, $rscript, $arg_hash);
	runProgram("mpiexec $executable", $arg_hash);
    }else{
	##doAnalysis($executable, $rscript, $arg_hash);
	runProgram($executable, $arg_hash);
    }
}
#subroutine to compare all files in sourcedir with originals and move
sub CompareThenMove {
    my $sourcedir = @_[0];
    my $prefix = "old_";
# define commands for Linux/UNIX
    my $diffcmd = "diff -s";
    my $movecmd = "mv -f";
    my $delcmd="rm -r";
# define command for Windows
    if($^O eq "MSWin32") {$diffcmd = "fc"; $slash="\\"; $movecmd = "move /y"; $delcmd="rmdir /s /q";}

    if (-e $sourcedir) {
      if (-e "$prefix$sourcedir") { # compare with sourcedir
	opendir(SOURCE, $sourcedir) or die "can't open $sourcedir folder: $!";
	while ( defined (my $file = readdir SOURCE) ) {
	  next if (($file =~ /^\.\.?$/) || ($file eq "logfile.txt"));     # skip . and .. and logfile
	  system("$diffcmd $sourcedir$slash$file $prefix$sourcedir$slash$file");
	  if($^O eq "MSWin32"){system("pause")}  # for checking comparisons in Windows
	} #compare
	closedir(SOURCE);

	my $olddir = "bk_$sourcedir";
	if (-e $olddir){system("$delcmd $olddir");} #delete old results directory (necessary for next command)
	system("$movecmd $prefix$sourcedir $olddir"); #preserve old results by renaming directory
      }

      system("$movecmd $sourcedir $prefix$sourcedir"); #rename results dir

    }
    else {print "$sourcedir does not exist\n"}
}

