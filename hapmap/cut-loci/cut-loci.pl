#!/usr/bin/perl -w
#
# Chops a big HAPMIXMAP file according to a small sample file.
#
# FIXME: Implement offsetting.

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

my $hm = Hapmix::HapmixData->new($data_base_name);
my $small = Hapmix::Genotypes->new($small_file);

# Test loci writing
# my $loci = Hapmix::Loci->new("chr7_phased_loci.txt");
# my @loci_names = $small->get_loci();
# $loci->write_file_from_loci_ids("all_loci.txt", $loci->get_loci_names());

# Columns to appear in the output file. Represented as column IDs.
my @loci_small = $small->get_loci();

# Write chosen columns (by IDs) to a file.
$hm->write_by_loci_ids($save_as, @loci_small);
