#!/usr/bin/perl -w
#
# Chops a big HAPMIXMAP file according to a small sample file.

use strict;

use Getopt::Long;
use Hapmix::HapmixData;
use Hapmix::Loci;

my $usage = 0;
my $data_base_name = '';
my $small_file = '';
my $help = 0;
my $offset = 0;
my $save_as = '';

GetOptions(
    "help!"  => \$usage,
    "data=s" => \$data_base_name,
    "small-file=s" => \$small_file,
    "offset=s" => \$offset,
    "save-as=s" => \$save_as,
);

$data_base_name or warn("Please specify data.\n");
$small_file or warn("Please specify small data file.\n");
$save_as or warn("Please specify --save-as basename.\n");

if (!$data_base_name or !$small_file or $help or !$save_as) {
    $usage = 1;
}

if ($usage) {
    print "\n";
    print "Usage: $0 [ options ]\n";
    print "   --help\n";
    print "   --data <basename>\n";
    print "   --small-file <filename>\n";
    print "   --offset <integer>\n";
    print "   --save-as <basename>\n";
    print "\n";
    exit(1);
}

# 1. Find the first and last column in sample, remember loci ids
# 2. Find loci in the locus file
# 3. Calculate the offset loci
# 4. Chop off column from the big file and rows from the loci file.
# 5. Write everything to new files.

my $hm = Hapmix::HapmixData->new($data_base_name);
my $small = Hapmix::Genotypes->new($small_file);

# print $small->get_headers() . "\n";
# print join(" ", $small->get_loci()) . "\n";
# my $first = $small->get_first_header();
# my $last = $small->get_last_header();

# my $loci = Hapmix::Loci->new("${data_base_name}_loci.txt");

# print $loci->{HEADER_LINE};
# print join("\n", $loci->get_lines_by_names(@list)) . "\n";

# print join("\n", $loci->get_lines_by_names($small->get_loci())) . "\n";
# print $small->{INDIVS} . "\n";
# print $small->get_columns($small->get_loci());

# Columns to appear in the out file. Represented as column IDs.
my @loci = $small->get_loci();

# Write chosen columns (by IDs) to a file.
$hm->write_genotypes_by_loci_ids("test.txt", @loci);
