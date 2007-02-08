#!/usr/bin/perl
# vim:set sw=4 softtabstop=4 ts=8 expandtab:

use strict;
use Getopt::Long;
use HapMix::Common;
use HapMix::Project;

my $project_name;
my $population;
my $usage = 0;
my $ok = 1;

GetOptions(
    "project=s" => \$project_name,
    "help!" => \$usage,
);

sub print_usage {
    print "\n";
    print "This script pre-processes the data.";
    print "\n";
#     print "Usage:\n";
#     print "$0 --name <project name> --population <pop>\n";
#     print "\n";
}

if (not $project_name) {
    print "Please provide a project name.\n";
    $ok = 0;
}

if (not $ok or $usage) {
    print_usage();
    exit(1);
}

my $project = HapMix::Project->new($project_name);

# * Find out which chromosome it is by reading the loci identifiers.
# * Order the loci.
# * Download data from hapmap.org.
# * Cut out the corresponding loci with given offset (may accept a parameter)
# * Write the genotypes and locus file for hapmixmap.
# * Write the ouctome file.

$project->find_chromosomes();
