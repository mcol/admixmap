 ** Instructions for running HAPMIXMAP with the example data set **

[information about HAPMIXMAP goes here]

[information about the data goes here]

[instructions for extracting the compressed file go here]

Running the program
-------------------------------------------------------
To run the program , you will need:
a genotypesfile
a locusfile
a text file containing a list of options

To run a case-control analysis you will also need:
a case-control genotypes file
an outcome variable file
a file with initial values of the allele frequencies
a file with the initial values of the parameters of the allele frequency prior
a file with initial values of the arrival rates

These last three are produced in a 'training' run, using only the HapMap data.

Examples of all these files are provided so you can try using the program right away.

A perl script is provided as a convenient way to run the program. To use this you must have perl installed.
If you do not have perl, you can still specify program options by writing them in a text file.

If you are using Linux/Unix, perl should already be installed. If you are using Windows, go to 
http://aspn.activestate.com/ASPN/Downloads/ActivePerl/ to get ActivePerl.

To run the script type 'perl hapmixmap.pl' into a console.

Program options are stored in the script using a 'hash'. Entries look like this:
 samples => 250,
 paramfile => "ParameterSamples.txt"

The options are then written to file and passed to the program.

Modify the options by editing the entries in the hash. Refer to the manual for details on program options.

You can specify the more commonly-changed option on the command-line. Type 'perl hapmixmap.pl --help' for a full list.
For example, you can specify the number of iterations and the length of the burnin by typing:
perl hapmixmap --samples=250 --burnin=50

If the HAPMIXMAP executable is not in the same directory as the script, you will need to change the line
my $executable = "./hapmixmap";

Make sure that the names of data file are correct.

To do an initial training run with only the HapMap data use the --init option.
To do further training runs, using the final parameter values from the previous run as initial values, omit --init.
This is not necessary for the example data as initial value files are already provided but if you wish to try it, make sure you keep a copy of the initial value files provided.


To run a case-control analysis, use the --cc option. This will tell the program to read the case-control genotypes file and the outcome variable file, use the saved values of the parameters as initial values and implement a score test (see manual for details of this).
The final parameter values will be saved in case you want to do further runs with the same data. You can change the names of these files by editing the appropriate option names.


Preparing Data
----------------------------------------------------------
The data files provided are already formatted for HAPMIXMAP but if you want to prepare your own data, this is how:

These instructions assume all typed loci lie on a single chromosome. Tools for preparing genome-wide data are in development.

Prepare a genotypes file and an outcomevarfile. See the manual for details on file formats. Ensure the genotype codings are consistent with those in the HapMap. 


Next obtain the HapMap data files for the appropriate chromosome from the following urls, replacing N in chrN with the number of the chromosome:
http://www.hapmap.org/downloads/phasing/2006-07_phaseII/phased/genotypes_chrN_CEU_r21_nr_fwd_legend.txt.gz
http://www.hapmap.org/downloads/phasing/2006-07_phaseII/phased/genotypes_chrN_CEU_r21_nr_fwd_phased.gz
http://www.hapmap.org/downloads/phasing/2006-07_phaseII/phased/genotypes_chrN_CEU_r21_nr_fwd_sample.txt.gz

Save these files in the same directory as chrN_legend.txt, chrN_phased.txt and chrN_sample.txt, respectively, again replacing N with the number of the chromosome.

A program called FormatHapMapData is provided to turn these files into HAPMIXMAP formatted data files. First you will need to compile this using a C++ compiler, any one will do. Run the program with the following command-line options:
-c<chromosome number>
-p<path to directory with the three files listed above>
-l<path and name of locusfile to output>
-g<path and name of genotypesfile to output>

use the -h option for a full list of user options. 

A perl script is provided to automate this process. You may want to edit some of the settings, such as the chromosome number, filenames, and the C++ compiler. To run this simply type 'perl preparedata.pl' in a console.
