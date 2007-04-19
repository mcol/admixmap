#!/usr/bin/perl -w
# 
# Reads phased Hapmix format file and writes two files.
# One file contains haploid data from most of the individuals.
# The CC (case-control) data is written as diploid data.

use strict;
use Getopt::Long;
use List::Util 'shuffle';

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

# Get two lines of an individual by individual no
sub get_indiv_lines {
	my $indiv_no = shift;
	my $lines = shift;
	my $idx = $indiv_no * 2;
	return ($lines[$idx], $lines[$idx + 1]);
}

# Arrays for individuals
my (@indivs_cc, @indivs_train);
@indivs_cc = @indiv_idx[0 .. ($no_masked - 1)];
@indivs_train = @indiv_idx[$no_masked .. $#indiv_idx];

sub compose_diploid(\@\@$) {
	# my $lines_ref = @_;
	# my @lines = @$lines_ref;
	my ($cols1_ref, $cols2_ref, $limit) = @_;
	my @outlist;
	my @cols1 = @$cols1_ref;
	my @cols2 = @$cols2_ref;
	$outlist[0] = shift(@cols1);
	shift(@cols2);
	foreach my $i (0 .. ($limit - 1)) {
		$outlist[$i + 1] = '"' . $cols1[$i] . "," . $cols2[$i] . '"';
	};
	return join("\t", @outlist) . "\n";
}

sub write_idx_to_file(\$\$$\@\@$) {
	my ($filename_ref, $header_ref, $format, $indices_ref, $wlines_ref, $limit_loci) = @_;
	my $filename = $$filename_ref;
	my $header = $$header_ref;
	my @indices = @$indices_ref;
	my @wlines = @$wlines_ref;
	my $limit = 0;
	# Valid formats
	my %formats = (
		'diploid' => 1,
		'haploid' => 1);
	if (!($formats{$format})) {
		die("Format ($format) should be either diploid or haploid.");
	}
	open(OUT_FILE, ">$filename");
	my @firstline = split(/\s+/, $header);
	if ($format eq "diploid") {
		$firstline[0] = "Individ";
	}
	$limit = $limit_loci ? $limit_loci : ($#firstline - 1);
    # If the limit is less than the number of loci, use the number of
    # loci.
    if ($limit_loci and $limit_loci < $#firstline) {
        $limit =  $limit_loci;
    } else {
        $limit =  ($#firstline - 1);
    }
	print OUT_FILE join("\t", @firstline[0 .. $limit]) . "\n";
	foreach my $idx (@indices) {
		my @ind_lines = get_indiv_lines($idx, @wlines);
		# print substr($ind_lines[0], 0, 72) . "\n";
		# print substr($ind_lines[1], 0, 72) . "\n";
		@cols1 = split(/\s+/, $ind_lines[0]);
		@cols2 = split(/\s+/, $ind_lines[1]);
		# Check for label correctness, every two rows:
		# LABEL
		# LABEL_2
		if (!($cols2[0] eq ($cols1[0] . "_2"))) {
			die("Wrong gamete labels: $cols1[0], $cols2[0].");
		}
		if ($format eq 'haploid') {
			print OUT_FILE join(" ", @cols1[0 .. $limit]) . "\n";
			print OUT_FILE join(" ", @cols2[0 .. $limit]) . "\n";
		} elsif ($format eq 'diploid') {
			print OUT_FILE compose_diploid(@cols1, @cols2, $limit);
			# print compose_diploid(@cols1, @cols2);
		} else {
			die("This shouldn't happen.");
		}
	}
	close(OUT_FILE);
}

sub rewrite_loci($$$) {
	my ($in, $out, $limit) = @_;
	open(IN_LOCUSFILE, "<$in")
		or die ("Can't open the input locus file '$in' for reading.");
	open(OUT_LOCUSFILE, ">$out")
		or die ("Can't open the output locus file '$out' for writing.");
	my $counter = 0;
	while (<IN_LOCUSFILE>) {
		if ($limit && $counter <= $limit) {
			print OUT_LOCUSFILE $_;
		}
		$counter++;
	}
	close(IN_LOCUSFILE);
	close(OUT_LOCUSFILE);
}

write_idx_to_file($ccFile,      $header, "diploid", @indivs_cc,    @lines, $limit_loci);
write_idx_to_file($haploidFile, $header, "haploid", @indivs_train, @lines, $limit_loci);
rewrite_loci($in_locus_file, $out_locus_file, $limit_loci);

