# CXX=pathCC
DFLAGS=

all:
	make -C .. all

#convertdata creates ADMIXMAP data files from HapMap data
convertdata: ../src/tools/convertdata.o
	$(CXX) -o convertdata ../src/tools/convertdata.o

#convertphaseddata creates ADMIXMAP data files from phased HapMap data
convertphaseddata: ../src/tools/convertphaseddata.o
	$(CXX) -o convertphaseddata ../src/tools/convertphaseddata.o

#sortSNPs sorts the SNPs in the HapMap data files by position (there are a number of loci out of sequence)
#run this before convertdata
sortSNPs: ../src/tools/sortSNPs.o
	$(CXX) -o sortSNPs ../src/tools/sortSNPs.o -lgsl -lgslcblas

#combi combines a case/control genotypes file with HapMap genotypes, setting genotypes at untyped loci to missing
#can also be used to set a subset of the HapMap genotypes to missing at a given set of loci
combi:
	$(CXX) $(DFLAGS) -ocombi makeIlluminaGenotypes.cc

clean:
	rm -f convertdata sort combi
