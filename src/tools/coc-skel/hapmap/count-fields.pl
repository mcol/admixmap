#!/usr/bin/perl -w

use strict;

# Counts fields in a whitespace-separated file.

my @els;

foreach my $line (<>) {
	# @els = split(/\t/, $line);
	@els = split(/\s+/, $line);
	print "Elements: " . scalar(@els) . "\n";
}
