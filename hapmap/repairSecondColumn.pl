#!/usr/bin/perl -w
# 
# author: Maciej Blizinski
#
# This script repairs the second column in data from the www.hapmap.org
# project. It reads from the standard input and writes to the standard
# output.
#
# Typical usage:
#
# ./repairSecondColumn.pl < infile > outfile

use strict;

my @cols;
my $line;
my @gts;

sub unique_gts {
	my @gts = @_;
	my @letters;
	my $two_letters;
	foreach $two_letters (@gts) {
		# Skip missing values
		if ($two_letters eq "NN") {
			next;
		}
		push(@letters, substr($two_letters, 0, 1));
		push(@letters, substr($two_letters, 1, 1));
	}
	# Unique values of @letters
	my %hash = map { $_ => 1 } @letters;
	return sort keys %hash;
}

my $i = 0;

# Read from the standard input.
while(<>) {
	$line = $_;
	# First line (headers) should be just printed as it is
	if ($i == 0) {
		print $line;
		$i++;
		next;
	};
	@cols = split(/\s+/, $_);
	my @indicated = split("/", $cols[1]);
	# Change only when there are more than 2 indicated genotypes,
	# for example G/C/T instead of actual G/C.
	if (@indicated > 2) {
		# Genotypes are from 11th column on
		$cols[1] = join("/", unique_gts(@cols[11 .. $#cols]));
	}
	# Print out the ready line to the standard output.
	print join(" ", @cols). "\n";
	$i++;
}
