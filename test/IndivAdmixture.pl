#!/usr/bin/perl -w


my $DEBUG = 1; # zero gives less output

# Change this to the location of the admixmap executable
my $executable = '../admixmap/admixmap';
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

    logfile                      => "IndResults/logfile.txt",
    paramfile                    => "IndResults/paramfile.txt",
    indadmixturefile             => "IndResults/indadmixture.txt"
};

doAnalysis($executable,$arg_hash);
&RenameDir("IndResults", "IndResults2Way");

$arg_hash->{initalpha0} = "1,1,1";
$arg_hash->{initalpha1} = "1,1,1";
doAnalysis($executable,$arg_hash);
&RenameDir("IndResults", "IndResults3Way");

$arg_hash->{fixedallelefreqs} = 1;
doAnalysis($executable,$arg_hash);
&RenameDir("IndResults", "IndResults3WayFixed");


sub dumpArgs
{
    my ($hash,$path) = @_;
    open(FILE,">$path") or die($!);
    foreach my $key (keys %$hash){
	print FILE $key ."\t" . $hash->{$key} ."\n";
    }
    close(FILE);
}

sub doAnalysis
{
    my ($prog,$args) = @_;
    my $command = $prog.getArguments($args);
    unless (-e "$results"){
	system("mkdir IndResults");
    }
    dumpArgs($args,"IndResults/args.txt");
    print $command,"\n" if $DEBUG;
    system($command);
    print "Starting R script to process output\n";
    system('RCMD BATCH --quiet --no-save --no-restore IndivAdmixture.R IndResults/Rlog.txt');
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
