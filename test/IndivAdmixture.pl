#!/usr/bin/perl -w

# Change this to the location of the admixmap executable
my $executable = './admixmap';
my $results = "IndResults";

# $arg_hash is a hash of parameters passed to
# the executable as arguments.
#
# keys (left-hand side) are parameter names
# values (right-hand side) are parameter values
my $arg_hash = 
{
    burnin   => 500, 
    samples  => 4500,
    every    => 10,
    locusfile                    => "IndData/loci.txt",
    genotypesfile                => "IndData/genotypes.txt",
    priorallelefreqfile          => "IndData/priorallelefreqs3way.txt",
    randommatingmodel            => 1,
    globalrho                    => 0,
    initalpha0                   => "1,1,1", # parameter vectors for Dirichlet prior on admixture 
    initalpha1                   => "1,1,0",
    chib                         => 1,
    #thermo                       => 1,
    #numannealedruns             => 100,

    resultsdir => "IndResults",
    logfile                      => "logfile.txt",
    paramfile                    => "paramfile.txt",
    indadmixturefile             => "indadmixture.txt"
};

$arg_hash->{resultsdir} = "IndResults2Way";
doAnalysis($executable,$arg_hash, "IndResults2Way");

$arg_hash->{initalpha0} = "1,1,1";
$arg_hash->{initalpha1} = "1,1,1";
$arg_hash->{resultsdir} = "IndResults3Way";
#doAnalysis($executable,$arg_hash, "IndResults3Way");

$arg_hash->{fixedallelefreqs} = 1;
$arg_hash->{resultsdir} = "IndResults3WayFixed";
#doAnalysis($executable,$arg_hash, "IndResults3WayFixed");


sub doAnalysis
{
    my ($prog,$args, $resultsdir) = @_;
    my $command = $prog.getArguments($args);
    unless (-e "$resultsdir"){
	system("mkdir $resultsdir");
    }
    system($command);
    print "Starting R script to process output\n";
     system("R --quiet --no-save --no-restore <AdmixmapOutput.R >results/Rlog.txt RESULTSDIR=$resultsdir");
    print "R script completed\n\n";
}

sub getArguments
{
    my $hash = $_[0];
    my $arg = '';
    foreach my $key (keys %$hash){
	$arg .= ' --'. $key .'='. $hash->{$key};
    }
    return $arg;
}

sub RenameDir{
    my($source, $target) = @_;
    if ($^O eq "MSWin32") {
	if (-e $target) {system("rmdir /s /q $target")}
	system("move /y $source $target");
    }
    else {
	if (-e $target) {system("rm -r $source")}
   	system("mv -f  $source $target");
    }
}
