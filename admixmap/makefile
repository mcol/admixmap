# makefile for compiling admixmap
CC 			= gcc
CXX 			= g++

MATRIXLIBDIR 		= ../gslwrap-ch/.libs/
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
CPPFLAGS		= -W -Wall -march=i686 $(IFLAGS) $(DFLAGS)

all:			admixmap .depend

again: clean admixmap

clean:          
		- rm -f del
		rm -f *.o
		rm -f .depend


veryclean:		
		- rm -f del
		rm -f *.o
		rm -f .depend
		rm -f ../test/admixmap
		rm -Rf $(DISTDIR)
		rm -f $(DISTDIR).tar
		rm -f $(DISTDIR).tar.gz
		rm -Rf $(WINDISTDIR)
		rm -f $(WINDISTDIR).zip

.depend: *.cc
	$(CXX) $(CPPFLAGS) -MM *.cc >> .depend

infiles		= admixmap.o Latent.o Regression.o AlleleFreqs.o AdmixOptions.o VectorLoop.o $(MATRIXLIBDIR)libmat.a DARS.o TuneRW.o gaussian_d.o Genome.o Chromosome.o CompositeLocus.o IndividualCollection.o Individual.o MetropolisHastings.o HaploidTransitionMatrix.o HMM.o IndAdmixOutputter.o AlleleFreqOutputter.o functions.o chib.o StratificationTest.o ScoreTests.o DispersionTest.o LogWriter.o InputData.o StringSplitter.o StringConvertor.o

admixmap:	$(infiles)
		$(CXX) $(CPPFLAGS) $(LFLAGS) -o ../test/admixmap $(infiles) $(LIBS)

check:		admixmap
		perl ../test/batchtest.pl

#Compile and make a tar.gz archive with data, docs and test scripts (Linux)
dist:		admixmap
		rm -Rf $(DISTDIR)
		mkdir $(DISTDIR)
		cp ../test/admixmap $(DISTDIR)
		cp README.txt $(DISTDIR)
		cp COPYING.txt $(DISTDIR)
		cp ../test/testArguments.txt $(DISTDIR)
		cp AdmixmapOutput.R $(DISTDIR)
		cp IndivAdmixture.R $(DISTDIR)
		cp admixmap.pl $(DISTDIR)
		cp batchtest.pl $(DISTDIR)
		cp IndivAdmixture.pl $(DISTDIR)
		mkdir $(DISTDIR)/data
		mkdir $(DISTDIR)/IndData
		mkdir $(DISTDIR)/results
		cp ../test/data $(DISTDIR)
		cp ../test/IndData $(DISTDIR)
		tar cf $(DISTDIR).tar $(DISTDIR)
		gzip --best $(DISTDIR).tar

#Compile  for Linux and repackage distribution with new exec
linuxdist:      admixmap
		cp ../test/admixmap $(DISTDIR)
		rm -f ../test/ $(DISTDIR).tar.gz
		tar cf $(DISTDIR).tar $(DISTDIR)
		gzip --best $(DISTDIR).tar

#Compile for Windows and prepare to repackage
#Zipping cannot be done in Windows from command line unless a special zip app installed (eg Winzip, AlZip etc)
windist:	admixmap
		cp ../test/admixmap.exe $(DISTDIR)
		rm -f ../test/ $(DISTDIR).zip

#command to zip directory - do it manually for now
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
include .depend
