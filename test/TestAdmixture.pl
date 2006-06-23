##
# Purpose: Given a genotypesfile, runs a single individual analysis in ADMIXMAP for each individual, 
# for each of two models, one with admixture, one without
#

print "script began: ";
my $starttime = scalar(localtime());
print $starttime;
print "\n";

# path to admixmap exec
my $execpath = "~"; 
if($^O eq "MSWin32"){
  $execpath =  "c:/cvs";
}
my $executable = "$execpath/genepi/test/admixmap";
my $genotypesfile = "data/genotypes.txt";
my $reportfile = "reportfile.txt";

#options for program
my $arg_hash = 
{
    burnin   => 5,
    samples  => 25   ,
    every    => 1,
    randommatingmodel            => 1,
    globalrho                    => 0,

 #chib => 1,
 thermo => 1,
 #numannealedruns => 100,
    #globalsumintensitiesprior    => '3,0.5'
    resultsdir => 'results',
 displaylevel => 1,
    locusfile                    => 'data/loci.txt',
    priorallelefreqfile          => 'data/priorallelefreqs.txt',
    logfile                      => 'log.txt',
    #ergodicaveragefile           => 'ergodicaverages.txt',
    indadmixturefile             => 'individualVarSamples.txt'
};
#open genotypesfile for input
open(GENOFILE, $genotypesfile) or die("Could not open genotypes file");

#specify file for individual genotypes
my $indivgenofile = "data/indivgenotypes.txt";
$arg_hash->{genotypesfile} = $indivgenofile;

#open report file for output and write header
print "writing report to $reportfile";
open(REPORTFILE, ">$reportfile") or die("Could not open report file");
print REPORTFILE "Individual\tNoAdmixtureModel\tAdmixtureModel\n";

#read header of genotypesfile
my $locusnames= <GENOFILE>;

my $index = 1;
##foreach my $line(<GENOFILE>){
while($index < 3){
$line = <GENOFILE>;
#open indiv genotypes file
  open(INDIVGENOFILE, ">$indivgenofile") or die("could not open indiv genotype file");
#write header
  print INDIVGENOFILE "$locusnames\n";
#write individual's genotypes
  print INDIVGENOFILE "$line";
  close(INDIVGENOFILE);

  #write line number
  print REPORTFILE "$index\t";

##Model with no admixture
  $arg_hash->{admixtureprior}         = "1,0,0,";
  $arg_hash->{admixtureprior1}        = "1,0,0,";
  doAnalysis($executable,$arg_hash);
  writeMargLikelihood($arg_hash);

##Model with admixture
  $arg_hash->{admixtureprior}         = "1,1,1,";
  $arg_hash->{admixtureprior1}        = "1,1,1,";
  doAnalysis($executable,$arg_hash);
  writeMargLikelihood($arg_hash);

  print REPORTFILE "\n";
  $index = $index +1;
}

close GENOFILE;
close REPORTFILE;
##move report into resultsdir
#system("mv $reportfile $arg_hash->{resultsdir}/$reportfile");
print "script ended: ";
my $endtime = scalar(localtime());
print $endtime;
print "\n";

sub writeMargLikelihood{
  my ($args) = @_;
  my $logfilepath = "$args->{resultsdir}/$args->{logfile}";
  open(LOGFILE, $logfilepath) or die("Could not open log file");

  foreach my $line (<LOGFILE>) {
    chomp($line);
    if($line =~ /ChibAlgorithm/ | $line =~ /LogEvidence/) {
      $line =~ s/[a-z]+//gi;
      print REPORTFILE "$line\t";
    }
  }
  close(LOGFILE);
}

sub getArguments {
    my $hash = $_[0];
    my $arg = '';
    foreach my $key (keys %$hash) {
	$arg .= ' --'. $key .'='. $hash->{$key};
    }
    return $arg;
}

sub doAnalysis {
    my ($prog, $args) = @_;
    my $command = $prog.getArguments($args);

    $ENV{'RESULTSDIR'} = $args->{resultsdir};
    #print "Results will be written to subdirectory $ENV{'RESULTSDIR'}";
    system($command);
    #print "Starting R script to process output\n";
    #print("R CMD BATCH --quiet --no-save --no-restore \
     #      ~/genepi/test/AdmixmapOutput.R $args->{resultsdir}/Rlog.txt\n");
    #system("R CMD BATCH --quiet --no-save --no-restore \
     #      ~/genepi/test/AdmixmapOutput.R $args->{resultsdir}/Rlog.txt");
    #print("Rcmd completed\n");
}
