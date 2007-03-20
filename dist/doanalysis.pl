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
    print "\nResults will be written to subdirectory $args->{resultsdir}\n";

    if($^O eq "MSWin32") {
	$prog =~ s/\//\\/g;
	if( !($prog =~ m/\.exe/) ) {
	    $prog = "$prog.exe";##check: should run ok without .exe extension
	}
    }
    if( !(-e $prog) ) {die "$prog not found\n"}; 

    runProgram($prog, $args);
    removeTrailingCommaFromTxtFiles($args->{resultsdir});
    runRScript($rscript, $args);
}

sub doParallelAnalysis {
    my ($prog, $rscript, $args) = @_;
    if( !(-e $prog) ) {die "$prog not found\n"}; 
    $ENV{'RESULTSDIR'} = $args->{resultsdir};
    runProgram("mpiexec $prog", $args);
    removeTrailingCommaFromTxtFiles($args->{resultsdir});
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

    print "Starting R script to process output\n";
    $ENV{'RESULTSDIR'} = $args->{resultsdir};
    system("$rcmd $rcmdArgs $rscript $args->{resultsdir}/Rlog.txt\n");
    print "R script completed\n\n";
}

sub removeTrailingCommaFromTxtFiles {
    my $dirname = $_[0];
    my $content;
    opendir(DIR, $dirname) or die "can't opendir $dirname: $!";
    my $file;
    while (defined($file = readdir(DIR))) {
        # do something with "$dirname/$file"
        if ($file =~ /.txt$/) {
            my $file_name = "$dirname/$file";
            $content = readFileToString($file_name);
            writeStringToFile($file_name, removeTrailingComma($content));
        }
    }
    closedir(DIR);
}

sub readFileToString {
    my $file_name = $_[0];
    # print "readFileToString(): Reading file '$file_name'\n";
    open(IN_FILE, "<$file_name")
        or die ("Couldn't open '$file_name' for reading.");
    my @lines = <IN_FILE>;
    close(IN_FILE);
    my $s = "@lines";
    return $s;
}

sub writeStringToFile {
    my $file_name = $_[0];
    my $s = $_[1];
    open(OUT_FILE, "<$file_name")
        or die ("Couldn't open '$file_name' for reading.");
    print OUT_FILE $s;
    close(OUT_FILE);
}

sub removeTrailingComma {
    my $s = $_[0];
    $s =~ s/,\s*\)/)/mg;
    return $s;
}

return (1);
