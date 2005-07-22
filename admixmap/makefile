# makefile for compiling admixmap
CC 			= gcc
CXX 			= g++

MATRIXLIBDIR 		= ../gslwrap-ch/
LIBDIR			= ../gslwrap-ch/

DISTDIR			= ../test/admixmap-linux
WINDISTDIR		= ../test/admixmap-win

IFLAGS			= -I$(LIBDIR)

LFLAGS			= -L$(MATRIXLIBDIR) -L$(LIBDIR)
LIBS			= -lmat -lm -static -lgsl -lgslcblas

# Release version flags
 DFLAGS = -O3

# Debug version flags
# DFLAGS = -g

# Profiling version flags
# DFLAGS = -pg -O3

CXXFLAGS		= 
CPPFLAGS		= -W -Wall $(IFLAGS) $(DFLAGS)

all:			admixmap .depend

again: clean admixmap

clean:          
		- rm -f del
		rm -f *.o
		rm -f .depend

.depend: *.cc
	$(CXX) $(CPPFLAGS) -MM *.cc >> .depend

infiles		= admixmap.o Latent.o Regression.o AlleleFreqs.o AdmixOptions.o $(LIBDIR)libmat.a DARS.o TuneRW.o gaussian_d.o Genome.o Chromosome.o CompositeLocus.o IndividualCollection.o Individual.o MetropolisHastings.o HMM.o IndAdmixOutputter.o functions.o chib.o StratificationTest.o ScoreTests.o MisSpecAlleleFreqTest.o HWTest.o DispersionTest.o LogWriter.o InputData.o StringSplitter.o StringConvertor.o DirichletParamSampler.o HamiltonianMonteCarlo.o

admixmap:	$(infiles)
		$(CXX) $(CPPFLAGS) $(LFLAGS) -o ../test/admixmap $(infiles) $(LIBS)

check:		admixmap
		perl ../test/batchtest.pl

include: .depend

#Compile  for Linux and repackage distribution with new exec
linuxdist:      admixmap
		cp ../test/admixmap $(DISTDIR)
		rm -f ../test/$(DISTDIR).tar.gz
		tar -czf admixmap.tar.gz ../test/$(DISTDIR)

#Compile for Windows and prepare to repackage
#Zipping cannot be done in Windows from command line unless a special zip app installed (eg Winzip, AlZip etc)
windist:	admixmap
		cp ../test/admixmap.exe $(WINDISTDIR)
		rm -f ../test/$(WINDISTDIR).zip
## ALZip
#		alzip -a ../test/$(WINDISTDIR) ../test/admixmap.zip
## WinZip
#		??
## pkzip
# 		pkzip admixmap.zip ../test/$(WINDISTDIR)


