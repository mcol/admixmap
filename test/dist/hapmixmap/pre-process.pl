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

# 1. Validate the data syntax.
# 2. Find out which chromosome it is by reading the loci identifiers.
# 3. Order the loci.
# 4. Download data from hapmap.org.
# 5. Cut out the corresponding loci with given offset (may accept
#    a parameter)
# 6. Write the genotypes and locus file for hapmixmap.
# 7. Write the ouctome file.

$project->validate();

