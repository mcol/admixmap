#!/usr/bin/perl
use strict; 
use File::Path;

sub getArguments
{
    my $hash = $_[0];
#    my $arg = '';
#   foreach my $key (keys %$hash){
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
    my $command = $prog.getArguments($args);
    
    $ENV{'RESULTSDIR'} = $args->{resultsdir};
    print "\nResults will be written to subdirectory $ENV{'RESULTSDIR'}\n";
    system($command);
    #print "Starting R script to process output\n";
    #system("R CMD BATCH --quiet --no-save --no-restore ../test/AdmixmapOutput.R $args->{resultsdir}/Rlog.txt");
    #print "R script completed\n\n";
}

sub CompareThenMove {
    my ($sourcedir, $targetdir) = @_;
    my $prefix = "old_";
# define commands for different OS's
    my $diffcmd = "diff -s"; my $movecmd = "mv -f"; my $slash="/";my $delcmd="rm -R";
    if (-e $targetdir) { # compare with sourcedir
	opendir(SOURCE, $sourcedir) or die "can't open $sourcedir folder: $!";
	while ( defined (my $file = readdir SOURCE) ) {
	    next if (($file =~ /^\.\.?$/) || ($file eq "logfile.txt"));     # skip . and .. and logfile
	    system("$diffcmd $sourcedir$slash$file $targetdir$slash$prefix$file");
	    if($^O eq "MSWin32"){system("pause")}  # for checking comparisons 
	} #compare
	closedir(SOURCE);
    }
    if (-e $sourcedir) {
	if (-e $targetdir) {
	    my $olddir = "$prefix$targetdir";
	    if (-e $olddir){system("$delcmd $prefix$targetdir");} #delete old results directory (necessary for next command)
	    system("$movecmd $targetdir $prefix$targetdir"); #preserve old results by renaming directory
	    }  
	system("$movecmd $sourcedir $targetdir"); #rename results dir
	mkdir("$sourcedir");                      #we need results dir for next analysis
	opendir(TARGET, $targetdir) or die "can't open $targetdir folder: $!";
	while ( defined (my $file = readdir TARGET) ) { #for each results file
	    next if $file =~ /^\.\.?$/;     # skip . and ..
	    system("$movecmd $targetdir$slash$file $targetdir$slash$prefix$file ");} #prefix with "old_"
	closedir(TARGET);
    }
    else {print "$sourcedir does not exist\n"}
}

################### DO NOT EDIT ABOVE THIS LINE ########################
my $executable = '~/test/admixmap';

my $arg_hash = {
#data files
    genotypesfile                   => 'data/genotypes.txt',
    locusfile                       => 'data/loci.txt',
    #priorallelefreqfile             => 'data/allelefreqs.txt',
    #fixedallelefreqs                => 1,
#main options
    samples  => 250,
    burnin   => 50,
    every    => 1,
    numannealedruns => 0,
    displaylevel => 3, 
    hapmixmodel => 1,
    hapmixlambdaprior=>"1000,10,40",

    #initialhapmixlambda => 0.4,
#sampler settings for lambda
    rhosamplerparams=> "0.01,  0.0001,  1.0,  0.9,  20",

#output file options
    resultsdir=> 'results',
    paramfile => 'paramfile.txt',
    allelefreqprioroutputfile => 'freqparams.txt',
# ergodicaveragefile => 'ergodicaverages.txt',
    logfile                     => 'log.txt',
# optional tests
    residualallelicassocscorefile => 'residualLDscoretests.txt',
    #hwscoretestfile                   => 'HardyWeinbergtest.txt'
};

# model with 4 block states
$arg_hash->{populations}           = 4;
$arg_hash->{resultsdir}            = 'Results';  
doAnalysis($executable,$arg_hash);
#CompareThenMove("Results", "Results4");

############### DO NOT EDIT BELOW THIS LINE ############################
