#!/usr/bin/perl -w
use strict;
use File::Path;
use Cwd;

#my $dir = getcwd;
my @Panels = ("CEU", "YRI", "JPTCHB");
my $whichchr = "chr22_5kloci";
#my $whichchr = "chr22";

my $dataprefix = "data/$whichchr/hapmixmap";
my $mytempdir = "/exports/work/scratch/pmckeigu";
my $resultsdir = "";

my @taskfiles = ("compare_states_train_tasks.txt",
		 "compare_states_test_tasks,txt", 
                 "compare_priors_train_tasks.txt",
		 "compare_priors_test_tasks,txt"); 
                 
#foreach my $taskfile(@taskfiles) {
#my $taskfile = "compare_priors_train_tasks.txt";
my $taskfile = "test_taskfile.txt";

open(TASKS, $taskfile) or die ("could not open $taskfile\n");

foreach my $command (<TASKS>) { # loop over lines of task file
    chomp($command);  # remove the newline from $line
    my $configfile = $command;
    $configfile =~ s/hapmixmap //g;
    #print "config file is $configfile\n";
    
    open (IN, "+<$configfile"); # open the config file with read/write access
    my @contents = <IN>; # contents of file held in array 
    seek IN,0,0; # set filehandle position
    foreach my $option (@contents) { # ? loop over lines of file
	$option =~ s/=results/=$mytempdir/; # substitute tempdir
	if($option =~ m/resultsdir=$mytempdir/) {
	    $resultsdir = $option;
            $resultsdir =~ s/\n//;
	    $resultsdir =~ s/resultsdir=//;
	    #print "resultsdir is $resultsdir\n";
	    if( !(-e $resultsdir) ) { # create resultsdir if required
		mkpath $resultsdir or die "cannot make directory $resultsdir";
	    }
	    unlink("$resultsdir/*");
	}
	print IN $option;
    }
    close IN;

    my $command = "hapmixmap.pathscale $configfile";
    my $script = $configfile;
    $script =~ s/.txt/.sh/;
    open(SCRIPT, "> $script");  
    print(SCRIPT "$command\n");
    close(SCRIPT);
    my $args = "-cwd -e $resultsdir.err -o $resultsdir.out -l h_rt=00:00:30 -N hapmix -V $configfile.sh";
    print "qsub $args\n";
    system("qsub $args");
    close TASKS;
};

