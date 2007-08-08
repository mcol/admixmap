Included in this directory are tools for downloading and formatting HapMap data as well as case-control genotypes files, for use with HAPMIXMAP.

Formatting unphased HapMap data
---------------------------------
Note: This program is obsolete and unsupported. HAPMIXMAP is intended to be used with phased HapMap data. This section is for information only.

You will have to download the unphased data yourself. For each population, you will need the information file, pedinfo2sample_XXX.txt.gz (where XXX is a the population code), currently located at http://www.hapmap.org/downloads/samples_individuals/.
Then, for each chromosome, you will need the 'non-redundant' fwd strand genotypes file, genotypes_chrY_XXX_r22_nr.bZZ_fwd.txt.gz (where Y is the chromosome number and ZZ is the current build).
You will then need to run the perl script repairSecondColumn.pl to ensure the foramatting program reads the bases correctly. Finally run the formatting program, FUD. Type FUD -h for a list of options.

Formatting phased HapMap data
--------------------------------
You can choose to either download the HapMap data files yourself or use the perl script getdata.pl to do it for you. 
For each chromosome and each population, there are 3 files, *.legend.txt.gz, *.sample.txt.gz and *.phased.gz.

Having downloaded and unzipped the data files, run the formatting program FPHD to produce a correctly-formatted locusfile and genotypesfile. The program can also be run via the getdata.pl script.
You can also supply a case-control genotypes file in the format described below and the program will not only encode the genotypes to match the HapMap encoding but also produce reduced locusfile and genotypesfile, containing only loci in the region typed. This will reduce redundancy and memory requirements when running HAPMIXMAP. It can also subdivide long chromosomes where necessary.
Type FPHD -h for a full list of program options. A typical usage would be as follows:
./FPHD -pHapMapData/CEU/chr7 -lHapMapData/CEU/chr7loci.txt -gHapMapData/CEU/chr7genotypes.txt -imyData/rawgenotypes.txt -omyData/CaseControlgenotypes.txt -maxloci=50000

Specifying maxloci, causes the chromosome to be subdivided if there are more than 50000 loci. To specify the size of the overlap between segments and the size of the flanking region either side of the typed region, use the 'minoverlap' and 'flanksize' options.

A case-control genotypes file for input to the formatting program should be as follows. The format is the same as the normal ADMIXMAP/HAPMIXMAP genotypesfile (sex columns are not yet supported) but with the genotypes coded as pairs of alleles (preferably fwd-strand coding but this is not essential). Note that the program assumes bases are never paired with their complements (ie A:T or C:G). Haploid genotypes are not yet supported.

idno	rsxxxx	rsxxxx	rsxxxx ...	
12345	G:A	A:C	C:C    ...
67890	G:A	A:A	C:A    ...	
...

To use the script you will need the Perl modules Getopt::Long and LWP::Simple. Both are available from CPAN (www.cpan.org) and usually part of the perl distribution. For a full list of options type 'perl getdata.pl -h' . The two required arguments are a chromosome number (-c) and a population code (-p). The script can be used to do either or both of of the downloading (-d) and formatting (-f).
A typical usage for both downloading and formatting might be as follows:
perl getdata.pl -c=7 -p=CEU -d -f -locusfile=HapMapData/CEU/chr7loci.txt -genotypesfile=HapMapData/CEU/chr7genotypes.txt -g=myData/rawgenotypes.txt -ccgfile=myData/CaseControlgenotypes.txt




