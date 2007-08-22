HAPMIXMAP running instructions.

1. Introduction
-------------------------------
HAPMIXMAP is a program for modelling extended haplotypes in genetic association studies. It is mainly intended to model HapMap haplotypes using tag SNP genotype data. It is released under a GPL (GNU public licence). See the file COPYING for details.

2. Requirements
-------------------------------
hapmixmap requires R to be installed on your computer. R is a free statistical analysis program that can be downloaded from http://www.r-project.org. 

You will need to add the path for R to your environment variable PATH. 

In Windows XP, go to Start-Control Panel-System. On the Advanced tab, click on Environment Variables. Select the PATH variable and click on Edit. Add the path of the R bin directory to the end of the path variable. This is usually something like ";C:\Program Files\R\bin".

In Linux/Unix use the command  "export PATH=$PATH:<path to R bin directory>". For convenience, this should be in the script that loads on login. 

Use of the Perl script requires Perl to be installed. It can be downloaded from http://www.perl.com/download.csp. 

3. Contents
------------
Distributed with the hapmixmap executable are the following:
  options.txt: a template options file with all possible options listed
  AdmixmapOutput.R: an R script to analyze the output. (Yes, it is the same as the one that goes with ADMIXMAP)
  README.txt:    this file.
  COPYING:       a copy of the GNU Public Licence (GPL).
  INSTALL:       instructions on installing from source
  data:          a directory with tutorial data
  tutorial.pl:   a perl script for the tutorial
  doanalysis.pl: a script required by the tutorial script
  training-initial.conf, training-resume.conf, testing.conf: 
                 options files for the tutorial
  FPHD:          a program to format HapMap data and case-control genotypes files
  getdata.pl:    a perl script to download data from HapMap and run the FPHD program  

4. Installing and Running
-------------------------------
To install from source, please refer to the INSTALL file.

Copy the package to a suitable location. Unzip/unpack the package. Navigate to the hapmixmap directory.
In Windows, unzip the package hapmixmap-x.xx.zip (in, say C:\hapmixmap). In a command window (Start->Run->cmd), navigate to the hapmixmap directory (cd c:\hapmixmap)

The recommended way to run the program is through a Perl script. There is an example for the tutorial included with this distribution.

Alternatively, the options may be entered on the command line (not recommended) or supplied in a text file. A template options text file,
options.txt, is supplied with this distribution.
To run using  options in a text file, navigate to the hapmixmap directory and type
./hapmixmap (or just hapmixmap in Windows) followed by the name of the options file.
Note: Do not use the name 'args.txt' for an options file as this is written to by the program.

Consult the documentation for details of the options.
Consult the tutorial documentation for more details on running the program.

You will know when the program has finished running when you see the word
Finished
followed by a line of *s

5. Formatting phased HapMap data
--------------------------------
You can choose to either download the HapMap data files yourself or use the perl script getdata.pl to do it for you. 
For each chromosome and each population, there are 3 files, *.legend.txt.gz, *.sample.txt.gz and *.phased.gz.

Having downloaded and unzipped the data files, run the formatting program FPHD to produce a correctly-formatted locusfile and genotypesfile. The program can also be run via the getdata.pl script.
You can also supply a  testgenotypes file in the format described below and the program will not only encode the genotypes to match the HapMap encoding but also produce reduced locusfile and genotypesfile, containing only loci in the region typed. This will reduce redundancy and memory requirements when running HAPMIXMAP. It can also subdivide long chromosomes where necessary.
Type FPHD -h for a full list of program options. A typical usage would be as follows:
./FPHD -pHapMapData/CEU/chr7 -lHapMapData/CEU/chr7loci.txt -gHapMapData/CEU/chr7genotypes.txt -imyData/rawgenotypes.txt -omyData/testgenotypes.txt -maxloci=50000

Specifying maxloci, causes the chromosome to be subdivided if there are more than 50000 loci. To specify the size of the overlap between segments and the size of the flanking region either side of the typed region, use the 'minoverlap' and 'flanksize' options.

A  testgenotypes file for input to the formatting program should be as follows. The format is the same as the normal ADMIXMAP/HAPMIXMAP genotypesfile (sex columns are not yet supported) but with the genotypes coded as pairs of alleles (preferably fwd-strand coding but this is not essential). Note that the program assumes bases are never paired with their complements (ie A:T or C:G). Haploid genotypes are not yet supported.

idno	rsxxxx	rsxxxx	rsxxxx ...	
12345	G:A	A:C	C:C    ...
67890	G:A	A:A	C:A    ...	
...

To use the script you will need the Perl modules Getopt::Long and LWP::Simple. Both are available from CPAN (www.cpan.org) and usually part of the perl distribution. For a full list of options type 'perl getdata.pl -h' . The two required arguments are a chromosome number (-c) and a population code (-p). The script can be used to do either or both of of the downloading (-d) and formatting (-f).
A typical usage for both downloading and formatting might be as follows:
perl getdata.pl -c=7 -p=CEU -d -f -locusfile=HapMapData/CEU/chr7loci.txt -genotypesfile=HapMapData/CEU/chr7genotypes.txt -g=myData/rawgenotypes.txt -ccgfile=myData/CaseControlgenotypes.txt

Troubleshooting
-----------------------
If you think you have found a bug or encounter any problems running the program or if it crashes or exits with an error message, contact the program authors. 
Contact information is available in the manual. If you specify invalid options or there is a problem in one of the data files, you should see a message stating what is wrong. To help diagnose any problem, please supply either the logfile or screen output as well as a list of the options you specified. We may also need a copy of the data in order to replicate the problem. If a problem occurs while running the R script, please supply the file Rlog.txt so we can see where.