#!/usr/bin/perl -w
# Remove trailing comma from R array, like this:
# a = c(1, 2, 3, )
#              ^
# Or with a new line
# a = c(1, 2, 3,
#              ^
# )
#
# a = c(1, 2, 3)
#
# Usage:
# ./remove-trailing-comma.pl file-name.txt

use strict;

sub remove_trailing_comma {
    my $s = $_[0];
    $s =~ s/,\s*\)/)/m;
    return $s;
}

my $file_name = $ARGV[0];

# Read the file
open(IN_FILE, "<$file_name") or die ("Couldn't open '$file_name' for reading.");
my @lines = <IN_FILE>;
my $all_text = "@lines";
close(IN_FILE);

# Write back to the (same) file.
open(OUT_FILE, ">$file_name") or die ("Couldn't open '$file_name' for writing.");
print OUT_FILE remove_trailing_comma($all_text);
close(OUT_FILE);
