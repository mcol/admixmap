#!/usr/bin/perl
use strict; 
use File::Path;
use File::Copy;
use File::Spec;
use File::Find;

sub getArguments {
    my $hash = $_[0];
    my $filename = 'perlargs.txt';
    open(OPTIONFILE, ">$filename") or die ("Could not open args file");
    foreach my $key (keys %$hash){
      print OPTIONFILE $key . '=' . $hash->{$key} . "\n";
    }
    close OPTIONFILE;
    return " ".$filename;
}

sub findFile{
  ##append .exe extension in Windows
  my $filename = $_[0];
  if($^O eq "MSWin32") {
    if( !($_[0] =~ m/\.exe/) ) {
      $filename = "$_[0].exe";
    }
  }

  ##first check if file exists and is a file (not a directory)
  if( (-e $filename) && (-f $filename) ){
    return("./$filename");
  }else{
    ##skip if filename is a path (contains '/')
    if(! ($filename =~ m(/)m)){

      ##search PATH variable
      my $slash = '/';
      #  if($^O eq "MSWin32") {
      #    $slash = '\\';
      #  }
      my @PATH = File::Spec->path();
      foreach(@PATH){
	my $path = "$_$slash$filename";
	if( (-e $path) && (-f $path)){
	  ##print "found $path\n";
	  return($path);
	}
      }
    }
  }
  return ('');
}

sub doAnalysis {
    my ($prog, $rscript, $args) = @_;

    my $path = findFile($prog);
    if( !$path ) {
      die "$prog not found\n"
    }

    my $resultsdir = $args->{resultsdir};
    if(!$resultsdir){
      $resultsdir = 'results';
    }
    print "\nResults will be written to subdirectory $resultsdir\n";
    if(runProgram($path, $args) == 0){
      runRScript($rscript, $args);
    }else{
      die("$prog failed to complete");
    }
}

sub doParallelAnalysis {
    my ($prog, $rscript, $args) = @_;
    if( !(-e $prog) ) {die "$prog not found\n"}; 

    runProgram("mpiexec $prog", $args);
    runRScript($rscript, $args);
}

sub runProgram {
    my ($prog, $args) = @_;
    system $prog.getArguments($args);
}

sub runRScript {
    my ($rscript, $args) = @_;
    my $rcmd = "R CMD";
    if($^O eq "MSWin32") {
	$rcmd = "Rcmd";##check: R CMD works on Windows too
	$rscript =~ s/\//\\/g;
    }

    if( !(-e $rscript) ) {die "$rscript not found\n"}; 
    my $rcmdArgs = "BATCH --quiet --no-save --no-restore";  

    my $resultsdir = $args->{resultsdir};
    if($resultsdir eq ''){
      $resultsdir = 'results';
    }
    print "Starting R script to process output\n";
    $ENV{'RESULTSDIR'} = $resultsdir;
    system("$rcmd $rcmdArgs $rscript $args->{resultsdir}/Rlog.txt\n");
    print "R script completed\n\n";
}


return (1);
