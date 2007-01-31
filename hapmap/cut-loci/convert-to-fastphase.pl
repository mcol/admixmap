#!/usr/bin/perl -w
#
# Converts haploid hapmixmap data into fastphase format.

use strict;

use Getopt::Long;
use Hapmix::HapmixData;
use Hapmix::Loci;

my $usage = 0;
my $in_base_name;
my $out_file;
my $no_loci_count;

GetOptions(
    "help!"  => \$usage,
    "haploid=s" => \$in_base_name,
    "fastphase=s" => \$out_file,
    "no-loci-count!" => \$no_loci_count,
);

if ($usage) {
    print "\n";
    print "Usage: $0 [ options ]\n";
    print "   --haploid <basename>    Basename for data files.\n";
    print "   --fastphase <file name> Output fastPHASE file name.\n";
    print "   --no-loci-count         Don't include the line with loci count.\n";
    print "\n";
    print "Or:\n";
    print "   --help                  See this message.\n";
    print "\n";
    exit(1);
 }

$out_file or die("Please specify --fastphase <file name>.\n");
$in_base_name or die("Please specify --haploid <base name>.\n");

my $hm = Hapmix::HapmixData->new($in_base_name);
$hm->write_fastphase($out_file, $no_loci_count);
