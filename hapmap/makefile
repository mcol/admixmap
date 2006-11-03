CC=pathCC

all: convertdata

convertdata:
	$(CC) -oconvertdata convertdata.cc

sortSNPs:
	$(CC) -osort sortSNPs.cc -lgsl-lgslcblas

combi:
	$(CC) -ocombi makeIlluminaGenotypes.cc

clean:
	rm -f convertdata sort combi
