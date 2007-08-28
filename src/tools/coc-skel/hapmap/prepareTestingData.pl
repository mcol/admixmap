#!/usr/bin/perl -w
# 
# Reads phased Hapmix format file and writes two files.
# One file contains haploid data from most of the individuals.
# The CC (case-control) data is written as diploid data.

use strict;
use Getopt::Long;
use List::Util 'shuffle';
use Hapmix::Masking;

my $infile = 0;
my $usage = 0;
my $haploidFile = "";
my $ccFile = "";
my $pcmasked = 12;
my $limit_loci = 0;
my $in_locus_file = '';
my $out_locus_file = '';

GetOptions(
	"phased-file=s" => \$infile,
	"haploid-file:s" => \$haploidFile,
	"case-control-file=s" => \$ccFile,
	"percent-indivs=i" => \$pcmasked,
	"limit-loci=i" => \$limit_loci,
	"in-locus-file=s" => \$in_locus_file,
	"out-locus-file=s" => \$out_locus_file,
	"help!" => \$usage);

$infile or $usage = 1;

if (not $haploidFile) {
	warn("Please specify --haploid-file.");
	$usage = 1;
}
if (not $ccFile) {
	warn("Please specify --case-control-file.");
	$usage = 1;
}

if ($usage) {
	print "\n";
	print "Usage: $0 --phased-file <filename> [ options ]\n";
	print "Options:\n";
	print "   --haploid-file <filename>      \n";
	print "   --case-control-file <filename> \n";
	print "   --percent-indivs <integer>     \n";
	print "   --limit-loci <integer>         \n";
	print "   --out-locus-file <filename>    \n";
	print "\n";
	print "Or: $0 --help\n";
	print "\n";
	exit(1);
}

open(PHASEDFILE, "<$infile") or die ("Can't open $infile for reading");
my @lines = <PHASEDFILE>;
close(PHASEDFILE);

my $header = shift(@lines);
my (@cols1, @cols2);
if (scalar(@lines) % 2 != 0) {
	die("Number of data rows should be even.");
}
my $indivs = scalar(@lines) / 2;
# Create an array of shuffled individual indices
my @indiv_idx = shuffle(0 .. ($indivs - 1));

# Number of masked individuals
my $no_masked = int($indivs * $pcmasked / 100.0);
print "No. CC individuals: $no_masked\n";
if ($no_masked < 1) {
	die("No individuals masked. Please check your options.\n");
}

write_idx_to_file($ccFile,      $header, "diploid", @indivs_cc,    @lines, $limit_loci);
write_idx_to_file($haploidFile, $header, "haploid", @indivs_train, @lines, $limit_loci);
rewrite_loci($in_locus_file, $out_locus_file, $limit_loci);

