#!/usr/bin/perl -w


my $DEBUG = 1; # zero gives less output

# Change this to the location of the admixmap executable
my $executable = '.admixmap';
my $results = "IndResults";

# $arg_hash is a hash of parameters passed to
# the executable as arguments.
#
# keys (left-hand side) are parameter names
# values (right-hand side) are parameter values
my $arg_hash = 
{
    burnin   => 1000, 
    samples  => 21000,
    every    => 5,
    analysistypeindicator     => -1, # -1 = single individual 
    coutindicator   => 0,
    locusfile                    => "IndData/loci.txt",
    genotypesfile                => "IndData/genotypes.txt",
    priorallelefreqfile          => "IndData/priorallelefreqs3way.txt",
    randommatingmodel            => 1,
    globalrho                    => 0,
    sumintensities	       => 99, # prior on sum of intensities parameter rho is gamma with shape 
                                         # parameter alpha, scale parameter beta 
                                         # 98 specifies alpha = beta = 0 (flat prior on log rho)
                                         # 99 specifies alpha = 1, beta = 0 (flat prior on rho)
                                         # if <90, value specifies alpha (with beta = 1)
                                         
    truncationpoint              => 15, # upper truncation point (lower truncation point is 1)
    # fixedallelefreqs             => 1,
    initalpha0                   => "1,1,0", # parameter vectors for Dirichlet prior on admixture 
    initalpha1                   => "1,1,0",

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
doAnalysis($executable,$arg_hash, "IndResults3Way");

$arg_hash->{fixedallelefreqs} = 1;
$arg_hash->{resultsdir} = "IndResults3WayFixed";
doAnalysis($executable,$arg_hash, "IndResults3WayFixed");


sub doAnalysis
{
    my ($prog,$args, $resultsdir) = @_;
    my $command = $prog.getArguments($args);
    unless (-e "$resultsdir"){
	system("mkdir $resultsdir");
    }
    print $command if $DEBUG;
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
