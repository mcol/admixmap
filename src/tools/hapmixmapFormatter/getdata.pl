#! /usr/bin/perl
#
#  getdata.pl
#  Perl script to download data files from HapMap.
#  and run them through the FPHD formatter to prepare for use with HAPMIXMAP.
#  Supports all FPHD options, including formatting of case-control genotypes files.
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

my $Verbose = 0;
my $Quiet = 0; # Not the oppposite of Verbose (yet)!
my $CHR = 0;
my $POP = "CEU";
my $Download = 0;
my $Unzip = 0;
my $Format = 0;
my $genotypesfile = '';
my $locusfile = '';
my $ccgenotypefile = '';
my $output_ccgfile='';
my $flank = 0;
my $limitloci = 0;
my $missing_string = 'N';
my $gzip_path = "/usr/bin";
#my $gzip_path = "c:\\msys\\bin\\";#for Windows users who have msys
my $HapMapPrefix = '';
## where to put downloaded files
my $dataprefix = '';
my $formatter_exec = "FPHD";
my $HelpNeeded = 0;

if(length(@ARGV) == 0){
  $HelpNeeded = 1;
}
GetOptions(
"h"               => \$HelpNeeded,
"v"               => \$Verbose,
"q"               => \$Quiet,
"d!"              => \$Download,
"u!"              => \$Unzip,
"f!"              => \$Format,
"c=i"             => \$CHR,
"p=s"             => \$POP,
"g=s"             => \$ccgenotypefile,
"gzip-path=s"     => \$gzip_path,
"hapmap-prefix=s" => \$HapMapPrefix,
"formatter=s"     => \$formatter_exec,
"prefix=s"        => \$dataprefix,
"ccgfile=s"       => \$output_ccgfile,
"genotypesfile=s" => \$genotypesfile,
"locusfile=s"     => \$locusfile,
"flank=f"         => \$flank,
"loci=i"          => \$limitloci,
"missing=s"       => \$missing_string
);

print "Unprocessed by Getopt::Long\n" if $ARGV[0];
foreach (@ARGV) {
  print "$_\n";
}

if($HelpNeeded){
print "\nUsage: getdata.pl <args> [options]\n\n";
print "Required arguments:\n";
print "-c=<1...22>             set chromosome number\n";
print "-p=<CEU|YRI|JPT+CHB|AS> set population. CEU = European; YRI = African;\n";
print "                                        AS = JPT+CHB = Asian\n";
print "\nOptions:\n";
print "-h                      print this help message and exit\n";
print "-v                      be verbose\n";
print "-q                      suppress output\n";
print "-d                      download files from HapMap\n";
print "-u                      unzip downloaded files (with gzip)\n";
print "-f                      format files\n";
print "-hapmix-prefix=<>       set HapMap prefix\n";
print "-gzip-path=<>           set path to gzip\n";
print "-formatter=<>           set path to formatting program\n";
print "-genotypesfile=<>       name of formatted HapMap genotypes file to write\n";
print "-locusfile=<>           name of HAPMIXMAP locusfile to write\n";
print "-g=<>                   name of case-control genotypes file to format\n";
print "-ccgfile=<>             name of formatted case-control genotypes file to write\n";
print "-flank=F                size of flanking region in Kb\n";
print "-loci=n                 maximum number of loci per chromosome \n";
print "                        (valid only if no case-control file specified)\n";
print "-missing                character denoting a missing genotype\n";
exit;
}
if(!$dataprefix){
  $dataprefix = "HapMapData$slash$POP";
}
##check arguments
if($CHR < 1 || $CHR > 22){
  die "Invalid chromosome number: $CHR \n";
}
if(!($POP eq "CEU" || $POP eq "YRI" || $POP eq "JPY+CHB" || $POP eq "AS")){
  die "Invalid population code: $POP \n";
}
#substitute AS (Asian)
if($POP eq "AS"){$POP = "JPT+CHB";}

if(!$Quiet){
##write info
  print "Using Settings:\n";
  print "------------------------------------\n";
  print "Chromosome = $CHR\n";
  print "Population = $POP\n";
  if($Download){
    print "Downloading files to $dataprefix\n";
  }
  else {print "Not downloading\n";}
  if($Format){
    print "Formatting with $formatter_exec\n";
    if($genotypesfile){print "writing genotypesfile: $genotypesfile\n";}
    if($locusfile){print "writing locusfile: $locusfile\n";}
    if($ccgenotypefile){
      write "also formatting user genotypes file: $ccgenotypefile\n";
      if($output_ccgfile){print "writing ccgenotypesfile: $output_ccgfile\n";}
      if($flank){print "flanking region: ${flank}Kb\n";}
    }
  }
  print "------------------------------------\n\n";
}

########################
## Script Starts Here ##
########################

## 1: download files from HapMap and decompress
if ($Download){
  if(!$HapMapPrefix){
    $HapMapPrefix = "http://www.hapmap.org/downloads/phasing/2006-07_phaseII/phased/genotypes_chr${CHR}_${POP}_r21_nr_fwd";
  }
  my $LegendFile = "${HapMapPrefix}_legend.txt.gz";
  my $PhasedFile = "${HapMapPrefix}_phased.gz";
  my $SampleFile = "${HapMapPrefix}_sample.txt.gz";

  my @HapMapFiles = ($LegendFile, $PhasedFile, $SampleFile);

  if(!(-e $dataprefix)){
    if(!$Quiet){print "making $dataprefix\n";}
    mkpath($dataprefix);
  }
  my $outLegendFile = "$dataprefix${slash}chr${CHR}_legend.txt";
  my $outPhasedFile = "$dataprefix${slash}chr${CHR}_phased.txt";
  my $outSampleFile = "$dataprefix${slash}chr${CHR}_sample.txt";
  my @outFiles = ($outLegendFile, $outPhasedFile, $outSampleFile);

  if(!$Quiet){print "Downloading files from HapMap\n";}
  for( my $i = 0; $i < 3; ++$i){
    my $zipfile = "@outFiles[$i].gz";
    my $content = getstore( @HapMapFiles[$i], $zipfile)
      or die "Couldn't get @HapMapFiles[$i]";
    if ($Unzip){
      #decompress using gzip (if available)
      system("${gzip_path}${slash}gzip -d $zipfile");
    }
  }

  if(!$Quiet){print "Download complete\n";}
}
## 2: Format files
if( $Format){
#  if(!$Quiet){print "Formatting data\n";}

  my $cmd = "$formatter_exec -c$CHR -p$dataprefix";
  if($genotypesfile){$cmd = $cmd . " -g$genotypesfile";}
  if($locusfile){$cmd  = $cmd . " -l$locusfile";}
  if($ccgenotypefile){
    $cmd = $cmd . " -i$ccgenotypefile";
    if($output_ccgfile){$cmd = $cmd . " -o$output_ccgfile";}
    if($flank){$cmd = $cmd . " -f$flank";}
  }else{
    if($limitloci){$cmd = $cmd . " -n$limitloci";}
  }
  if($Verbose){$cmd = $cmd . " -v";}
  $cmd = $cmd . " -m$missing_string";

  #print "$cmd";
  my $status = system($cmd);
  if($status == 0){
 #  if(!$Quiet){print "Formatting complete\n"; }
  }  else{
   print "Formatter exited with status $status\n";
  }
}

#if(!$Quiet){
  print "\nscript complete\n\n"
#}
