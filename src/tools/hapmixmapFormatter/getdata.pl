#! /usr/bin/perl
#
#  getdata.pl
#  Perl script to download data files from HapMap.
#  and run them through the FPHD formatter to prepare for use with HAPMIXMAP.
#
#   This program is free software distributed WITHOUT ANY WARRANTY.
#   You can redistribute it and/or modify it under the terms of the GNU General Public License,
#   version 2 or later, as published by the Free Software Foundation.
#   See the file COPYING for details.
#
## Copyright(c) David O'Donnell 2007
###############################################################################################
use strict;
use File::Path;
use LWP::Simple;
use Getopt::Long;

#############
## Settings #
#############
my $slash = "/";
if($^O eq "MSWin32"){ $slash = "\\";}

my $beQuiet = 0;
my $CHR = 0;
my $Panel = "";
my $Download = 0;
my $Unzip = 0;
my $Format = 0;
my $gzip_path = "/usr/bin";
if($^O eq "MSWin32"){
  $gzip_path = "c:\\msys\\bin\\";#for Windows users who have MSys
}
my $HapMapURL = "http://www.hapmap.org/downloads/phasing/2006-07_phaseII/phased";
## where to put downloaded files
my $dataprefix = '';
my $formatter_exec = "FPHD";
my $formatter_args='';
my $HelpNeeded = 0;

if(length(@ARGV) == 0){
  $HelpNeeded = 1;
}
GetOptions(
"h"               => \$HelpNeeded,
"q"               => \$beQuiet,
"d!"              => \$Download,
"u!"              => \$Unzip,
"f!"              => \$Format,
"c=i"             => \$CHR,
"p=s"             => \$Panel,
"gzip-path=s"     => \$gzip_path,
"hapmap-url=s"    => \$HapMapURL,
"formatter=s"     => \$formatter_exec,
"prefix=s"        => \$dataprefix,
"args=s"          => \$formatter_args

);

print "Unprocessed by Getopt::Long\n" if $ARGV[0];
foreach (@ARGV) {
  print "$_\n";
}

if($HelpNeeded || !$CHR || !$Panel){
print "\nUsage: getdata.pl <args> [options]\n\n";
print "Use this script to download data files from HapMap, decompress them\n"; 
print " and run them through the FPHD formatter to prepare for use with HAPMIXMAP.\n\n";
print "Required arguments:\n";
print "-c=<1...22>             chromosome number\n";
print "-p=<CEU|YRI|JPT+CHB|AS> panel code. (AS = JPT+CHB)\n";
print "\nOptions:\n";
print "-h                      print this help message and exit\n";
print "-q                      suppress output\n";
print "-d                      download files from HapMap\n";
print "-u                      unzip downloaded files (with gzip)\n";
print "-f                      format files\n";
print "-hapmap-url=<>          HapMap url \n";
print "-gzip-path=<>           path to gzip\n";
print "-formatter=<>           path to formatting program\n";
print "-args=\"\"                arguments to formatting program\n";
exit;
}
my $dataprefix_arg = $dataprefix;
if(!$dataprefix){
  $dataprefix = "HapMapData$slash$Panel";
  ##prefix arg to formatter - always forward slash
  $dataprefix_arg = "HapMapData/$Panel";
}
##check arguments
if($CHR < 1 || $CHR > 22){
  die "Invalid chromosome number: $CHR \n";
}
if(!($Panel eq "CEU" || $Panel eq "YRI" || $Panel eq "JPT+CHB" || $Panel eq "AS")){
  die "Invalid population code: $Panel \n";
}
#substitute AS (Asian)
if($Panel eq "AS"){$Panel = "JPT+CHB";}

if(!$beQuiet){
##write info
  print "Using Settings:\n";
  print "------------------------------------\n";
  print "Chromosome = $CHR\n";
  print "Panel      = $Panel\n";
  print "Data directory = $dataprefix\n";
  if($Download){
    print "Downloading from $HapMapURL\n";
    if($Unzip){
      print "Decompressing with ${gzip_path}${slash}gzip\n";
    }else{
      print "Not decompressing\n";
      $Format = 0;
    }
  }
  else {print "Not downloading\n";}
  if($Format){
    print "Formatting with $formatter_exec\n";
  }else {print "Not formatting\n"}
  print "------------------------------------\n\n";
}

########################
## Script Starts Here ##
########################

## 1: download files from HapMap and decompress
if ($Download){
  my $HapMapPrefix = "$HapMapURL/genotypes_chr${CHR}_${Panel}_r21_nr_fwd";

  my $LegendFile = "${HapMapPrefix}_legend.txt.gz";
  my $PhasedFile = "${HapMapPrefix}_phased.gz";
  my $SampleFile = "${HapMapPrefix}_sample.txt.gz";
  my @HapMapFiles = ($LegendFile, $PhasedFile, $SampleFile);

  ##create data dir if it doesn't exist
  if(!(-e $dataprefix)){
    if(!$beQuiet){print "making $dataprefix\n";}
    mkpath($dataprefix);
  }
  my $outLegendFile = "$dataprefix${slash}chr${CHR}_legend.txt";
  my $outPhasedFile = "$dataprefix${slash}chr${CHR}_phased.txt";
  my $outSampleFile = "$dataprefix${slash}chr${CHR}_sample.txt";
  my @outFiles = ($outLegendFile, $outPhasedFile, $outSampleFile);

  if(!$beQuiet){print "Downloading files from HapMap\n";}
  for( my $i = 0; $i < 3; ++$i){
    my $zipfile = "@outFiles[$i].gz";
    my $content = getstore( @HapMapFiles[$i], $zipfile);
      ##getstore returns the HTTP response code
      ## 200 = success
      ## 4xx, 5xx, 6xx = error
    if($content != 200){
      die("Unable to retrieve data from HapMap website")
    }

    if ($Unzip){
      #decompress using gzip (if available)
      system("${gzip_path}${slash}gzip -d $zipfile");
    }
  }

  if(!$beQuiet){print "Download complete\n";}
}
## 2: Format files
if( $Format){
#  if(!$beQuiet){print "Formatting data\n";}

  my $cmd = "$formatter_exec -c=$CHR -p=$dataprefix_arg $formatter_args";

  my $status = system($cmd);
  print "\nFormatter exited with status $status\n";
}

#if(!$beQuiet){
  print "\nscript complete\n\n"
#}
