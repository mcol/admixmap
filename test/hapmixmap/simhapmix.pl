#!/usr/bin/perl
use strict; 
use File::Path;
use Getopt::Long;

my $parallel = '';
my $simulate = '';
my $samples  = 25;
my $burnin = 5;
my $every = 1;
my $executable = '';
my $rscript = "../../dist/AdmixmapOutput.R";

my $serial_executable = "$ENV{'HOME'}/usr/bin/hapmixmap";
my $parallel_executable = "$ENV{'HOME'}/usr/bin/hapmixmap-para";

##settings for R script
##TODO: allow these to be set on the command-line
my $datadir = "simdata";
my $NumIndivs = 60;
my $states = 8;
my $NumCases = 10;
my $NumControls = 10;
my $chrLength = 0.5;
my $ExtremeAlleleFreqs = 0;
my $FixedMixtureProps = 1;
my $FreqDispersionShape = 2;
my $FreqDispersionRate = 10;

GetOptions("exec=s"    => \$executable,
           "rscript=s" => \$rscript,
           "parallel"  => \$parallel, 
	   "simulate"  => \$simulate, 
	   "samples=i" => \$samples,
	   "burnin=i"  => \$burnin,
	   "every=i"   => \$every);

if($simulate){
##write settings file for R script
  open(RINPUTFILE, ">simHapMixSettings.R") or die("Could not open Rinputfile file");
  print RINPUTFILE "datadir<-\"$datadir\"\n";
  print RINPUTFILE "N<-$NumIndivs\n";
  print RINPUTFILE "K<-$states\n";
  print RINPUTFILE "NumCases<-$NumCases\n";
  print RINPUTFILE "NumControls<-$NumControls\n";
  print RINPUTFILE "chr.L<-$chrLength\n";
  print RINPUTFILE "numChr<-length(chr.L)\n";
  if($ExtremeAlleleFreqs){
    print RINPUTFILE "extreme.allele.freqs<-T\n";
  }
  print RINPUTFILE "freq.dispersion.prior.shape<-$FreqDispersionShape\n";
  print RINPUTFILE "freq.dispersion.prior.rate<-$FreqDispersionRate\n";

  if($FixedMixtureProps){
    print RINPUTFILE "fixed.mixture.props<-T";
  }else{
    print RINPUTFILE "fixed.mixture.props<-F";
  }
  close(RINPUTFILE);

##run R simulation script
    print "Running R script to simulate data\n";
    system("R CMD BATCH --vanilla simHapMix.R simRlog.txt");
    print "simulation complete\n";
}

if(!$executable){
  if($parallel){
    $executable = $parallel_executable;
  }else{
    $executable = $serial_executable;
  }
}
my $function_file = "../../dist/doanalysis.pl";#"$ENV{'HOME'}/genepi/trunk/dist/doanalysis.pl";

require $function_file or die("cannot find doanalysis.pl");

################### DO NOT EDIT ABOVE THIS LINE ########################


my $arg_hash = {
#data files
    locusfile                       => "$datadir/loci.txt",

    #genotypesfile                   => "$datadir/genotypes.txt", #diploid data
     genotypesfile                   => "$datadir/genotypes_haploid.txt",#haploid data

    #priorallelefreqfile             => "datadir/allelefreqs.txt',
    #fixedallelefreqs                => 1,

#main options
    resultsdir => 'results',
    displaylevel   => 3, 

    samples  => $samples,
    burnin   => $burnin,
    every    => $every,

    numannealedruns => 0,
    thermo => 0,
    checkdata=> 0,
##model
    states=>$states,
    hapmixmodel => 1,
    fixedmixtureprops => $FixedMixtureProps,

##priors
arrivalrateprior=>"400, 1, 10, 1",
allelefreqprecisionprior => "$FreqDispersionShape, $FreqDispersionRate, 1",
freqprecisionhiermodel => 1,

arrivalratesamplerparams => "0.5, 0.00001, 10, 0.9, 20",

#output files
    logfile                     => 'logfile.txt',
    paramfile               => 'paramfile.txt',
    freqprecisionfile => "allelefreqpriorsamples.txt",
    #regparamfile          => 'regparamfile.txt',
    ergodicaveragefile => 'ergodicaverage.txt',

#final values
finalvaluedir => "initialValues",

#posterior means
    arrivalrateposteriormeanfile => "lambdaPosteriorMeans.txt",
    allelefreqprecisionposteriormeanfile => "freqDispersionPosteriorMeans.txt"

#optional tests
    #residualldtest => 1,
    #allelicassociationtest  => 1
};

# Initial run 
#haploid data
$arg_hash->{resultsdir}            = 'Results_Haploid';
callDoAnalysis();

#diploid data
#$arg_hash->{resultsdir}            = 'ResultsDiploid';
#$arg_hash->{genotypesfile} = "data/genotypes.txt";
#doAnalysis($executable,$arg_hash);
#system("cp Results/initialallelefreqs.txt data");

# rerun with final values of previous run as intial values of this
#system("cp $arg_hash->{resultsdir}/initialallelefreqs.txt data");
$arg_hash->{resultsdir}            = 'ResultsRerun';
$arg_hash->{initialvaluedir}="initialValues";
#$arg_hash->{fixedallelefreqs} = 0;
#delete $arg_hash->{priorallelefreqfile};
#doAnalysis($executable,$arg_hash);

sub callDoAnalysis {
    if($parallel){
	doParallelAnalysis($executable, $rscript, $arg_hash);
    }else{
	doAnalysis($executable, $rscript, $arg_hash);
    }
}
