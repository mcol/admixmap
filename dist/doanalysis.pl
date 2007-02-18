#!/usr/bin/perl
use strict; 
use File::Path;
use File::Copy;

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

sub doAnalysis {
    my ($prog, $rscript, $args) = @_;
    $ENV{'RESULTSDIR'} = $args->{resultsdir};
    print "\nResults will be written to subdirectory $ENV{'RESULTSDIR'}\n";
    my $rcmd = "R CMD";
    if($^O eq "MSWin32") {
	$rcmd = "Rcmd";
	$prog =~ s/\//\\/;
	$rscript =~ s/\//\\/;
    }
    if( !(-e $prog) ) {die "$prog not found\n"}; 
    if( !(-e $rscript) ) {die "$rscript not found\n"}; 
    my $rcmdArgs = "BATCH --quiet --no-save --no-restore";  
    system $prog.getArguments($args);
    print "Starting R script to process output\n";
    system("$rcmd $rcmdArgs $rscript $args->{resultsdir}/Rlog.txt\n");
    print "R script completed\n\n";
}
