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

Troubleshooting
-----------------------
If you think you have found a bug or encounter any problems running the program or if it crashes or exits with an error message, contact the program authors. 
Contact information is available in the manual. If you specify invalid options or there is a problem in one of the data files, you should see a message stating what is wrong. To help diagnose any problem, please supply either the logfile or screen output as well as a list of the options you specified. We may also need a copy of the data in order to replicate the problem. If a problem occurs while running the R script, please supply the file Rlog.txt so we can see where.