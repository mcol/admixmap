#!/usr/bin/perl
# script to test most ADMIXMAP options and compare results with previous results


use Getopt::Std;


use constant IS_MS		=> $^O eq "MSWin32"			;
use constant DIR_SEP		=> IS_MS ? "\\" : "/"			;
use constant DIFF_CMD		=> IS_MS ? "fc" : "diff"		;
use constant DIFFS_CMD		=> IS_MS ? "fc" : "diff -s"		;
use constant DEF_RDIR		=> "results"				;
use constant DEF_EXEC		=> "./admixmap"				;
use constant DEF_DDIR		=> "../../dist/admixmap/tutorial/data"	;
use constant PROG		=> ($0 =~ m/.*?([^\/]+)$/)[ 0 ]		;
use constant ARGS_FN		=> 'perlargs.txt'			;
use constant TEST_CATTLE	=> 0					;


$PROG = PROG; # variable is apparently easier to embed in strings than constant


#---------------------------------------------------------------------------
# Command-line options and usage errors:
#---------------------------------------------------------------------------

sub usage
    {
    print STDERR "\nUsage: $PROG <options>\n\n" .
		"Options:\n" .
		"	-h			 [print this message]\n" .
		"	-x <admixmap-executable> [defaults to " . DEF_EXEC . "]\n" .
		"	-d <data-dir>		 [defaults to " . DEF_DDIR . "]\n" .
		"	-r <results-dir>	 [defaults to " . DEF_RDIR . "]\n" .
		"	-o <output-file>	 [for admixmap's stdout, will be moved to results-dir]\n" .
		"	-E			 [admixmap's stderr is also in <output-file>]\n" .
		"	-C			 [compare admixmap's output against past runs\n" .
		"	-s			 [report results matching the last run (not just differences)]\n" .
		"	-e			 [prints commands to stdout before executing]\n" .
		"	-a			 [abort after error]\n" .
	    "The -o option likely will not work on MS-Windows.\n\n";
    exit 1
    }


my %args;
getopts( "hx:r:d:o:ECsea", \%args ) or usage;


usage if defined $args{"h"};


my $executable = defined $args{"x"} ? $args{"x"} : DEF_EXEC;
my $resultsdir = defined $args{"r"} ? $args{"r"} : DEF_RDIR;
my $datadir    = defined $args{"d"} ? $args{"d"} : DEF_DDIR;
my $outfile    = $args{"o"} if defined $args{"o"};
my $redir_err  = defined $args{"E"};
my $cmp_out    = defined $args{"C"};
my $diff_cmd   = defined $args{"s"} ? DIFFS_CMD : DIFF_CMD;
my $echo_cmds  = defined $args{"e"};
my $abort_err  = defined $args{"a"};

print	"\n\n**** BATCH TEST ****\n" .
	"\tOS: $^O\n" .
	"\tExecutable: $executable\n" .
	"\tData directory: $datadir\n" .
	"\tResults directory: $resultsdir\n";
print "\tAdmixmap's stdout" . ($redir_err?"/stderr":"") . " will be redirected into $outfile\n" if defined $outfile;
print "\n";



#---------------------------------------------------------------------------
# Execute a command, echoing if necessary.
#---------------------------------------------------------------------------

sub exec_cmd
    {
    my $command = shift( @_ );
    print PROG . ": execute: $command\n" if $echo_cmds;
    return system( $command );
    }



#---------------------------------------------------------------------------
# doAnalysis()
#---------------------------------------------------------------------------

sub doAnalysis
    {
    my ( $exec, $args, $run_desc ) = @_;

    print "\n==== Test run: $run_desc\n";

    # One way to test if results directory exists:
    if ( opendir $rdir,$resultsdir )
	{
	closedir $rdir;
	}
    else
	{
	print "$PROG: making directory \"$resultsdir\"\n";
	mkdir $resultsdir or die "$PROG: can't make directory $resultsdir: $!";
	}

    my $command = $exec . getArguments($args);

    if ( defined $outfile )
	{
	$command = $command . " >>$outfile";
	$command = $command . " 2>&1" if $redir_err;
	}

    my $status = exec_cmd( "$command" );
    if ( $status == -1 )
	{
	print STDERR PROG . ": error executing: $command : $!\n";
	}
    elsif ( $status != 0 )
	{
	print STDERR "Test failed: $run_desc\n";
	print STDERR PROG . ": $command exited with non-zero status: $status\n";
	}
    else
	{
	# Strip quotes that the "old" admixmap incorrectly leaves on locus
	# labels, making false differences with the new version's output.  This
	# would presumably be much better done within this perl script itself.
	exec_cmd( "sed -i -r 's/(^[^ \\t]+[ \\t]+[0-9]+[ \\t]+[0-9.]+[ \\t]+)\"([^\"]+)\"\$/\\1\\2/' \"$resultsdir"
		    . DIR_SEP . "LocusTable.txt\"" );
	}
    exec_cmd( "mv \"$outfile\" \"$resultsdir".DIR_SEP."$outfile\"" ) if ( defined $outfile );

    die if ( $abort_err && ($status != 0) );

    return ($status == 0);
    }


sub getArguments
{
    my $hash = $_[0];
#    my $arg = '';
#    foreach my $key (keys %$hash){
#	$arg .= ' --'. $key .'='. $hash->{$key};
#    }
#    return $arg;

    my $filename = $resultsdir . DIR_SEP . ARGS_FN;

    open(OPTIONFILE, ">$filename") or die ("Could not open args file \"$filename\"");
    foreach my $key (keys %$hash){
      print OPTIONFILE $key . '=' . $hash->{$key} . "\n";
    }
    close OPTIONFILE;

    if ( $echo_cmds )
	{
	print "\n\n\n======== OPTIONS FILE: ========\n";
	system( "cat " . " $filename" );
	print "===============================\n";
	}

    return " ".$filename;
}



#-----------------------------------------------------------------------------
# Execute a diff command
#-----------------------------------------------------------------------------

sub diff_files
    {
    my ( $file1, $file2 ) = @_;

    my $command = $diff_cmd . " $file1 $file2";

    exec_cmd( $command );

    system( "pause" ) if ($^O eq "MSWin32");	# for checking comparisons
						# (why no pause on non-MS?)
    }



#-----------------------------------------------------------------------------
# SUBROUTINE TO COMPARE ALL FILES IN sourcedir WITH ORIGINALS AND MOVE
#-----------------------------------------------------------------------------

sub CompareThenMove {
    my ($sourcedir, $targetdir) = @_;
    my $prefix = "old_";
# define commands for different OS's
    if($^O eq "MSWin32") {$movecmd = "move /y"; $delcmd="rmdir /s /q";}
    else {$movecmd = "mv -f"; $delcmd="rm -r";}
    if (-e $targetdir) { # compare with sourcedir
	opendir(SOURCE, $sourcedir) or die "can't open $sourcedir folder: $!";
	while ( defined (my $file = readdir SOURCE) ) {
	    next if (($file =~ /^\.\.?$/) || ($file eq "logfile.txt"));		   # skip . and .. and logfile
	    next if ( (defined $outfile) && (! $cmp_out) && ($file eq $outfile) ); # Don't compare amm's output
	    diff_files( "$sourcedir".DIR_SEP."$file", "$targetdir".DIR_SEP."$prefix$file" );
	    if($^O eq "MSWin32"){system("pause")}  # for checking comparisons
	} #compare
	closedir(SOURCE);
    }
    else
	{
	print "$PROG: **** no prior results with which to compare\n";
	}

    if (-e $sourcedir) {
	if (-e $targetdir) {
	    my $olddir = "$prefix$targetdir";
	    if (-e $olddir){exec_cmd("$delcmd $prefix$targetdir");} #delete old results directory (necessary for next command)
	    exec_cmd("$movecmd $targetdir $prefix$targetdir"); #preserve old results by renaming directory
	    }
	exec_cmd("$movecmd $sourcedir $targetdir"); #rename results dir
	mkdir("$sourcedir");                      #we need results dir for next analysis
	opendir(TARGET, $targetdir) or die "can't open $targetdir folder: $!";
	while ( defined (my $file = readdir TARGET) ) { #for each results file
	    next if $file =~ /^\.\.?$/;     # skip . and ..
	    exec_cmd("$movecmd $targetdir".DIR_SEP."$file $targetdir".DIR_SEP."$prefix$file ");} #prefix with "old_"
	closedir(TARGET);
    }
    else {print STDERR PROG . ": source-directory '$sourcedir' does not exist\n"}
}



#---------------------------------------------------------------------------
# List of test-runs to be performed:
#---------------------------------------------------------------------------

# $arg_hash is a hash of parameters passed to
# the executable as arguments.
#
# keys (left-hand side) are parameter names
# values (right-hand side) are parameter values
my $arg_hash = {
    samples                    => 15,
    burnin                     => 5,
    every                      => 1,
    locusfile                  => "$datadir/loci.txt",
    genotypesfile              => "$datadir/genotypes.txt",
    outcomevarfile             => "$datadir/outcomevars.txt",
    covariatesfile             => "$datadir/covariates3std.txt",
    displaylevel             => 2,
    "use-pedigree-for-individual" => 1,
    "no-conjugate-update"	  => 1,

# output files
    logfile                    => 'logfile.txt',
    paramfile                  => 'param.txt',
    regparamfile               => 'regparam.txt',
    indadmixturefile           => 'indadmixture.txt',
    ergodicaveragefile         => 'ergodicaverage.txt',

# extra output files
    haplotypeassociationtest  => 1,
    allelicassociationtest    => 1,
    allelefreqoutputfile      => 'allelefreqoutput.txt'
};

# single population, thermodynamic
$arg_hash->{thermo} = 1;
$arg_hash->{numannealedruns} = 4;
$arg_hash->{populations} = 1;
if(doAnalysis( $executable, $arg_hash, 'single population, thermodynamic' )){
    &CompareThenMove( $resultsdir, $resultsdir . '0' );
}

# single population, reference prior on allele freqs, annealing
$arg_hash->{thermo} = 0;
$arg_hash->{indadmixhiermodel} = 1;
$arg_hash->{stratificationtest}  = 1;
if(doAnalysis( $executable, $arg_hash, 'single population, reference prior on allele freqs, annealing' )){
    &CompareThenMove( $resultsdir, $resultsdir . '1' );
}

# two populations, reference prior on allele freqs
$arg_hash->{numannealedruns} = 0;
$arg_hash->{populations}     = 2;
if(doAnalysis( $executable, $arg_hash, 'two populations, reference prior on allele freqs' )){
    &CompareThenMove( $resultsdir, $resultsdir . '2' );
}

# fixed allele freqs, individual sumintensities
# possible problem here - changing numannealedruns should not change output except for annealmon.txt
delete $arg_hash->{populations};
delete $arg_hash->{allelefreqoutputfile};
$arg_hash->{fixedallelefreqs} = 1;
$arg_hash->{priorallelefreqfile}  = "$datadir/priorallelefreqs.txt",
$arg_hash->{allelefreqtest}  = 1;
$arg_hash->{allelefreqtest2} = 1;
$arg_hash->{ancestryassociationtest} = 1;
$arg_hash->{affectedsonlytest}       = 1;
$arg_hash->{globalrho} = 0;
$arg_hash->{numannealedruns} = 5;
$arg_hash->{thermo} = 0;
#if(doAnalysis( $executable, $arg_hash, 'fixed allele freqs, individual sumintensities' )){
    #&CompareThenMove( $resultsdir, $resultsdir . '3' );
#}

# prior on allele freqs
$arg_hash->{testoneindiv} = 0;
$arg_hash->{numannealedruns} = 0;
$arg_hash->{thermo} = 0;
$arg_hash->{fixedallelefreqs} = 0;
$arg_hash->{globalrho}        = 1;
$arg_hash->{randommatingmodel} = 1;
$arg_hash->{allelefreqoutputfile}     = 'allelefreqoutput.txt';
delete $arg_hash->{allelefreqtest};
delete $arg_hash->{allelefreqtest2};
delete $arg_hash->{affectedsonlytest};
$arg_hash->{dispersiontest}  = 1;
$arg_hash->{outcomevarcols} = 2; # skin reflectance
if(doAnalysis($executable,$arg_hash, 'prior on allele freqs')){
    &CompareThenMove( $resultsdir, $resultsdir . '4' );
}

# dispersion model for allele freqs
delete $arg_hash->{priorallelefreqfile};
delete $arg_hash->{dispersiontestfile};
$arg_hash->{historicallelefreqfile} = "$datadir/priorallelefreqs.txt";
$arg_hash->{affectedsonlytest}       = 1;
$arg_hash->{fstoutput} = 1;
$arg_hash->{dispparamfile} = 'disppar.txt';
$arg_hash->{randommatingmodel} = 0;
$arg_hash->{outcomevarcols} = 1; # diabetes
if(doAnalysis( $executable, $arg_hash, 'dispersion model for allele freqs' )){
    &CompareThenMove( $resultsdir, $resultsdir . '5' );
}

# prior on allele freqs, testoneindiv, no regression, thermo
# possible problem here - changing numannealedruns should not change output except for annealmon.txt
delete $arg_hash->{historicallelefreqfile};
delete $arg_hash->{affectedsonlytest};
delete $arg_hash->{ancestryassociationtest};
delete $arg_hash->{allelicassociationtest};
delete $arg_hash->{haplotypeassociationtest};
delete $arg_hash->{fstoutput};
delete $arg_hash->{dispparamfile};
delete $arg_hash->{regparamfile};
delete $arg_hash->{outcomevarcols};
delete $arg_hash->{outcomevarfile};
delete $arg_hash->{covariatesfile};
$arg_hash->{testoneindiv} = 1;
$arg_hash->{priorallelefreqfile}  = "$datadir/priorallelefreqs.txt",
$arg_hash->{randommatingmodel} = 1;
$arg_hash->{globalrho} = 0;
$arg_hash->{testoneindiv} = 1;
$arg_hash->{numannealedruns} = 100;
$arg_hash->{thermo} = 1;
if(doAnalysis($executable,$arg_hash,'prior on allele freqs, testoneindiv, no regression, thermo')){
    &CompareThenMove( $resultsdir, $resultsdir . '6' );
}


# For the moment, removed the cattle dataset, as it causes an error in admixmap:
if ( TEST_CATTLE ) {
    # autosomal and X chromosome data
    my $arg_hash = {
	genotypesfile             => 'Cattledata/ngu_complete2.txt',
	locusfile                 => 'Cattledata/loci.txt',
	populations		  => 2,
	globalrho		  => 1,
	samples			  => 25,
	burnin			  => 5,
	every			  => 1,
	resultsdir		  => "$resultsdir",
	logfile                   => 'log.txt',
	paramfile		  => 'paramfile.txt',
	indadmixturefile	  => 'indadmixture.txt',
	ergodicaveragefile	  => 'ergodicaverages.txt',
	allelefreqoutputfile	  => 'allelefreqs.txt',
	#allelicassociationtest   => 1,
	#ancestryassociationtest  => 1,
	#affectedsonlytest        => 1,
	#haplotypeassociationtest => 1,
	#stratificationtest       => 1
    };

    if(doAnalysis( $executable, $arg_hash, 'autosomal and X chromosome data' )){
        &CompareThenMove( $resultsdir, "cattleresults");
    }
}


# Single individual
my $arg_hash = {
    burnin   => 10,
    samples  => 60,
    every    => 1,
    numannealedruns => 0,
    displaylevel   => 2,

    locusfile                    => "IndData/loci.txt",
    genotypesfile                => "IndData/genotypes.txt",
    priorallelefreqfile          => "IndData/priorallelefreqs3way.txt",
    randommatingmodel            => 1,
    globalrho                    => 0,
    fixedallelefreqs             => 1,
    admixtureprior               => "1,1,0",
    admixtureprior1              => "1,1,1",
    logfile                      => "logfile.txt",
    chib                         => 1,
    indadmixturefile             => "indadmixture.txt"
};
if(doAnalysis( $executable, $arg_hash, 'Single individual' )){
    &CompareThenMove( $resultsdir, "Indresults");
}

print "\n";
