
Installing Admixmap for Windows
-------------------------------
Open the file admixmap.zip and unzip the contents into a suitable directory.

At the moment Admixmap cannot cope with filenames or directories containing spaces, therefore the contents of the zip file must be extracted to a directory which does not contain any spaces. I would therefore recommend creating the directory C:\Admixmap and extracting the contents to it.

For Admixmap to run correctly the directory structure MUST be preserved. Do NOT move the executable admixmap.exe or copy it to another location.

Installing Admixmap for Linux
-------------------------------
Unpack the file admixmap.tar.gz in a suitable directory (tar -xzf admixmap.tar.gz).

At the moment Admixmap cannot cope with filenames or directories containing spaces, therefore the contents of the zip file must be extracted to a directory which does not contain any spaces.

For Admixmap to run correctly the directory structure MUST be preserved. Do NOT move the executable admixmap or copy it to another location.

Requirements
------------
Admixmap requires R to be installed on your computer as well as the boa package. R is a free statistical analysis program that can be downloaded from http://www.r-project.org. The boa (Bayesian Output Analysis) package of functions for analysis of Markov chain Monte Carlo output can be installed from R with the command, 'install.package("boa")' or downloaded from the R website. 

If you are using Windows, you will need to add the path for R to your environment variable PATH. For Windows XP, go to Start-Control Panel-System. On the Advanced tab, click on Environment Variables. Select the PATH variable and click on Edit. Add the path of the R bin directory to the end of the path variable. This is usually something like ";C:\Program Files\R\bin".

Use of the Perl script requires Perl to be installed. It can be downloaded from http://www.perl.com/download.csp. 

Contents
------------
Distributed with the Admixmap executable are the following:
data		 a folder with test data and input files.
inddata		 a folder with data on a single individual.
testArguments.txt a text file containing a list of options to test the programme using the test data.
admixmap.pl	 a Perl script to test the programme using the test data.
indadmixture.pl	 a Perl script to demonstrate using the programme for inference about a
                 single individual.
admixoutput.R	 an R script to analyze the output.
indadmixture.R	 an R script to analyze the output when using data on an individual.
README.txt	 this file.
COPYING		 a copy of the GNU Public Licence (GPL).

Running the Programme
----------------------
The recommended way to run the programme is through a Perl script. There are two included with this distribution.
The options may be changed by editing these scripts. Consult the documentation for details of the options.
Alternatively, the options may be entered on the command line (not recommended) or supplied in a text file. A sample options text file,
testArguments.txt, with the same options as in admixmap.pl, is supplied with this distribution.

The test dataset is a cross-sectional sample of Hispanic Americans in San Luis Valley, California, typed at 21 marker loci, with measurements of a binary variable (diabetes) and a continuous variable (skin pigmentation, measured as skin reflectance).

The individual data is of a single (unknown) individual with no phenotype data.

To run the Perl scripts, open a console window(Start -> Run ->cmd in Windows), navigate to the admixmap directory and type
perl admixmap.pl or
perl indadmixture.pl
Depending on how Perl is installed, it may also be possible to run the scripts by double-clicking on them. 

To run using  options in a text file, navigate to the admixmap directory and type
./admixmap followed by the name of the options file.
Note: Do not use the name 'args.txt' for an options file as this is written to by the programme.