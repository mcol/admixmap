Introduction
-----------------
ADMIXMAP is a general-purpose program for modelling admixture using marker genotypes and trait data.  

The current version is 3.7 and it is released under a GPL (GNU public licence). See the file COPYING for details.

Installing Admixmap for Windows
-------------------------------
Open the file admixmap.zip and unzip the contents into a suitable directory.

At the moment Admixmap cannot cope with filenames or directories containing spaces, therefore the contents of the zip file must be extracted to a directory which does not contain any spaces. We therefore recommend creating the directory C:\Admixmap and extracting the contents to it.

For Admixmap to run correctly the directory structure MUST be preserved. Do NOT move the executable admixmap.exe or copy it to another location.

Installing Admixmap for Linux
-------------------------------
Unpack the file admixmap.tar.gz in a suitable directory (tar -xzf admixmap.tar.gz). You may need to change the permissions in order to run the program. Use the command 'chmod a=x admixmap' to do this.

At the moment Admixmap cannot cope with filenames or directories containing spaces, therefore the contents of the zip file must be extracted to a directory which does not contain any spaces.

For Admixmap to run correctly the directory structure MUST be preserved. Do NOT move the executable admixmap or copy it to another location.

Requirements
------------
Admixmap requires R to be installed on your computer. R is a free statistical analysis program that can be downloaded from http://www.r-project.org. 

You will need to add the path for R to your environment variable PATH. 

In Windows XP, go to Start-Control Panel-System. On the Advanced tab, click on Environment Variables. Select the PATH variable and click on Edit. Add the path of the R bin directory to the end of the path variable. This is usually something like ";C:\Program Files\R\bin".

In Linux/Unix use the command  "export PATH=$PATH:<path to R bin directory>". For convenience, this should be in the script that loads on login. 

Use of the Perl script requires Perl to be installed. It can be downloaded from http://www.perl.com/download.csp. 

Contents
------------
Distributed with the admixmap executable are the following:
  testArguments.txt: a text file containing a list of options to test the programme using the test data.
  admixmap.pl: a sample perl script to use with the test data.
  AdmixmapOutput.R: an R script to analyze the output.
  README.txt: this file.
  COPYING: a copy of the GNU Public Licence (GPL).
  tutorial: a folder with test data and tutorial scripts.

Running the Program
----------------------
The recommended way to run the program is through a Perl script. There are two included with this distribution.
The options may be changed by editing these scripts. Consult the documentation for details of the options.
Alternatively, the options may be entered on the command line (not recommended) or supplied in a text file. A sample options text file,
testArguments.txt, with the same options as in admixmap.pl, is supplied with this distribution.

The test dataset is a cross-sectional sample of Hispanic Americans in San Luis Valley, California, typed at 21 marker loci, with measurements of a binary variable (diabetes) and a continuous variable (skin pigmentation, measured as skin reflectance).

To run the Perl scripts, open a console window(Start -> Run ->cmd in Windows), navigate to the admixmap directory and type
perl admixmap.pl
Depending on how Perl is installed, it may also be possible to run the scripts by double-clicking on them. 

To run using  options in a text file, navigate to the admixmap directory and type
./admixmap followed by the name of the options file.
Note: Do not use the name 'args.txt' for an options file as this is written to by the programme.

You will know when the program has finished running when you see the word
Finished
followed by a line of *s

Troubleshooting
-----------------------
If you think you have found a bug or encounter any problems running the program or if it crashes or exits with an error message, contact the program authors. 
Contact information is available in the manual. If you specify invalid options or there is a problem in one of the data files, you should see a message stating what is wrong. To help diagnose any problem, please supply either the logfile or screen output as well as a list of the options you specified. We may also need a copy of the data in order to replicate the problem. If a problem occurs while running the R script, please supply the file Rlog.txt so we can see where.