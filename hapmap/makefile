# CC=pathCC
DFLAGS=

all: convertdata

#convertdata creates ADMIXMAP data files from HapMap data
convertdata:
	$(CC) $(DFLAGS) -oconvertdata convertdata.cc

#convertphaseddata creates ADMIXMAP data files from phased HapMap data
convertphaseddata:
	$(CC) $(DFLAGS) -oconvertphaseddata convertphaseddata.cc

#sortSNPs sorts the SNPs in the HapMap data files by position (there are a number of loci out of sequence)
#run this before convertdata
sortSNPs:
	$(CC) $(DFLAGS) -osort sortSNPs.cc -lgsl-lgslcblas

#combi combines a case/control genotypes file with HapMap genotypes, setting genotypes at untyped loci to missing
#can also be used to set a subset of the HapMap genotypes to missing at a given set of loci
combi:
	$(CC) $(DFLAGS) -ocombi makeIlluminaGenotypes.cc

clean:
	rm -f convertdata sort combi
