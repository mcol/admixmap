#!/usr/bin/perl
# vim:set sw=4 softtabstop=4 ts=8 expandtab:

use strict;
use Getopt::Long;
use HapMix::Common;
use HapMix::Project;

my $project_name;
my $population;
my $source_file;
my $usage = 0;
my $ok = 1;

my $common_data = HapMix::Common->new();

GetOptions(
    "name=s" => \$project_name,
    "population=s" => \$population,
    "genotypes=s" => \$source_file,
    "help!" => \$usage,
);

sub print_usage {
    print "\n";
    print "Usage:\n";
    print "$0 --name <project name> --population <pop>\n";
    print "\n";
}

if (not $project_name) {
    print "Please provide a project name.\n";
    $ok = 0;
}

if (not $population) {
    print "Please give a population name.\n";
    $ok = 0;
}

if (not $common_data->{POPULATIONS}{$population}) {
    print "Population can be one of ('";
    print join("', '", sort keys %{$common_data->{POPULATIONS}});
    print "'), '$population' was provided.\n";
    $ok = 0;
}

if (not $ok or $usage) {
    print_usage();
    exit(1);
}

if (-d "$project_name") {
    print "Directory $project_name already exists. Exiting.\n";
    exit(1);
}

my $project = HapMix::Project->new($project_name);
$project->initialize($population);
if ($source_file) {
    # print "Source file: $source_file\n";
    $project->get_source($source_file);
} else {
    print "Please place the source data in " . $project->get_source_file_name() . "\n";
}

