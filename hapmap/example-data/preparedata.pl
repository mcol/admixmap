#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Path;
use Time::Local;

my $HapMapDataDir = "data";
my $UserDataDir = "data";
my $GenotypesFile = "CaseControlGenotypes.txt";
my $Chromosome = "7";
my $HapMapDataPrefix = "$HapMapDataDir/chr$Chromosome";

#files with the formatted HapMap data
my $HapMapGenotypesFilename = "HapMapChr7Genotypes.txt";
my $HapMapLocusFilename = "HapMapChr7Loci.txt";

#directory with formatting tools
my $PrepDir = "prep";

#C++ compiler to compile formatting program
my $CC="g++";
#compiler flags
my $CPPFLAGS="-O3";

################## Script Starts Here ##################################

#compile formatting code
# print "Compiling formatting tools...\n";
# system("$CC $CPPFLAGS -o$PrepDir/FormatHapMapData $PrepDir/FormatHapMapData.cc");
# print "Compilation complete\n\n";

#convert HapMap data to HAPMIXMAP format
print "Converting HapMap Data to HAPMIXMAP format...\n";
system("$PrepDir/FormatHapMapData -c$Chromosome -p$HapMapDataDir -l$HapMapDataDir/$HapMapLocusFilename -g$HapMapDataDir/$HapMapGenotypesFilename");
print "Conversion completed\n\n";


print "Data are now ready for HAPMIXMAP\n\n";
