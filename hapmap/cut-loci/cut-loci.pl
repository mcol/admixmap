#!/usr/bin/perl -w
#
# Chops a big HAPMIXMAP file according to a small sample file.
#
# Range and offsetting implemented. Offsetting is quite inefficient, but
# is acceptable for now. Cutting test data: chromosome 7 with
# helen_genotypes.txt takes 4 minutes on Celeron 1.5GHz.

use strict;

use Getopt::Long;
use Hapmix::HapmixData;
use Hapmix::Loci;

my $default_offset = 100000; # == 100kb, kilo-base

my $usage = 0;
my $data_base_name = '';
my $sample_file = '';
my $help = 0;
my $save_as = '';
my $range = 0;
my $offset = $default_offset;

GetOptions(
    "help!"  => \$usage,
    "data=s" => \$data_base_name,
    "sample-file=s" => \$sample_file,
    "offset=s" => \$offset,
    "save-as=s" => \$save_as,
    "range!" => \$range,
);

$data_base_name or warn("Please specify data.\n");
$sample_file or warn("Please specify sample data file.\n");
$save_as or warn("Please specify --save-as basename.\n");

if (!$data_base_name or !$sample_file or $help or !$save_as) {
    $usage = 1;
}

if ($usage) {
    print "\n";
    print "Usage: $0 [ options ]\n";
    print "   --help\n";
    print "   --data <basename>       Basename for data files.\n";
    print "                           _genotypes.txt will be appended automagically.\n";
    print "   --sample-file <filename> Sample genotypes file.\n";
    print "                           This should contain _genotypes.txt.\n";
    print "   --offset <integer>      Offset in bp, $default_offset by default.\n";
    print "   --save-as <basename>    Basename for saved files\n";
    print "   --range                 Instead of taking specified loci, take a whole\n";
    print "                           range between the first and the last locus\n";
    print "                           in the sample.\n";
    print "\n";
    exit(1);
}

my $hm = Hapmix::HapmixData->new($data_base_name);
my $sample = Hapmix::Genotypes->new($sample_file);

# Test loci writing
# my $loci = Hapmix::Loci->new("chr7_phased_loci.txt");
# my @loci_names = $sample->get_loci();
# $loci->write_file_from_loci_ids("all_loci.txt", $loci->get_loci_names());

# Columns to appear in the output file. Represented as column IDs.
my @loci_sample = $sample->get_loci();

# If requested, consider given locus list as a range.
if ($range) {
    @loci_sample = $hm->range_by_ids(@loci_sample);
}

# If requested, apply offsetting
if ($offset) {
    @loci_sample = $hm->offset($offset, @loci_sample);
}

# Write chosen columns (by IDs) to a file.
$hm->write_by_loci_ids($save_as, @loci_sample);
