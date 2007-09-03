Included in this directory are tools for downloading and formatting HapMap data as well as test (eg case-control) genotypes files, for use with HAPMIXMAP.

Formatting phased HapMap data
--------------------------------
You can choose to either download the HapMap data files yourself or use the perl script getdata.pl to do it for you. 
For each chromosome and each population, there are 3 files, *.legend.txt.gz, *.sample.txt.gz and *.phased.gz.

Having downloaded and unzipped the data files, run the formatting program FPHD to produce a correctly-formatted locusfile and genotypesfile. The program can also be run via the getdata.pl script.
You can also supply a case-control genotypes file in the format described below and the program will not only encode the genotypes to match the HapMap encoding but also produce reduced locusfile and genotypesfile, containing only loci in the region typed. This will reduce redundancy and memory requirements when running HAPMIXMAP. It can also subdivide long chromosomes where necessary.
Type FPHD -h for a full list of program options. A typical usage would be as follows:
./FPHD -p=HapMapData/CEU/chr7 -l=HapMapData/CEU/chr7loci.txt -c=Chr7 -g=HapMapData/CEU/chr7genotypes.txt -i=myData/rawgenotypes.txt -o=myData/CaseControlgenotypes.txt -maxloci=50000

Specifying maxloci, causes the chromosome to be subdivided if there are more than 50000 loci. To specify the size of the overlap between segments and the size of the flanking region either side of the typed region, use the 'minoverlap' and 'flanksize' options.

A test genotypes file for input to the formatting program should be as follows. The format is the same as the normal ADMIXMAP/HAPMIXMAP genotypesfile (sex columns are not yet supported) but with the genotypes coded as pairs of alleles (preferably fwd-strand coding but this is not essential). Note that the program assumes bases are never paired with their complements (ie A:T or C:G). Haploid genotypes are not yet supported.

idno	rsxxxx	rsxxxx	rsxxxx ...	
12345	G:A	A:C	C:C    ...
67890	G:A	A:A	C:A    ...	
...

To use the script you will need the Perl modules Getopt::Long and LWP::Simple. Both are available from CPAN (www.cpan.org) and usually part of the perl distribution. For a full list of options type 'perl getdata.pl -h' . The two required arguments are a chromosome number (-c) and a population code (-p). The script can be used to do either or both of of the downloading (-d) and formatting (-f).
A typical usage for both downloading and formatting might be as follows:
perl getdata.pl -c=7 -p=CEU -d -u -f -args="-l=HapMapData/CEU/chr7loci.txt -g=HapMapData/CEU/chr7genotypes.txt -i=myData/rawgenotypes.txt -o=myData/CaseControlgenotypes.txt -maxloci=50000"



