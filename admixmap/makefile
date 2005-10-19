# makefile for compiling admixmap
CC 			= gcc
CXX 			= g++

DISTDIR			= ../test/admixmap-linux
WINDISTDIR		= ../test/admixmap-win

LIBS			= -lm -static -lgsl -lgslcblas

# processor architecture flags
I686FLAGS    = -march=i686 
P4FLAGS      = -march=pentium4 -mfpmath=sse
X86_64FLAGS  = -march=athlon64

# flags for release version, debug version or profiling version
# release 
DFLAGS = -O3
# debug  
# DFLAGS = -g
# profiling 
# DFLAGS = -pg -O3

CXXFLAGS		= 
CPPFLAGS		= -W -Wall $(I686FLAGS) $(DFLAGS)
# CPPFLAGS		= -W -Wall $(P4FLAGS) $(DFLAGS)

all:		rand admixmap .depend

again: clean all

clean:          
		- rm -f del
		rm -f *.o
		rm -f .depend

.depend: *.cc
	$(CXX) $(CPPFLAGS) -MM *.cc >> .depend

infiles		= admixmap.o Latent.o Regression.o AlleleFreqs.o AdmixOptions.o rand_serial.o AdaptiveRejection.o StepSizeTuner.o Gaussian.o Genome.o Chromosome.o CompositeLocus.o IndividualCollection.o Individual.o GaussianProposalMH.o HMM.o IndAdmixOutputter.o functions.o chib.o StratificationTest.o ScoreTests.o MisSpecAlleleFreqTest.o DispersionTest.o HWTest.o LogWriter.o InputData.o StringSplitter.o StringConvertor.o DirichletParamSampler.o HamiltonianMonteCarlo.o DataMatrix.o MuSampler.o DispersionSampler.o

rand:		
		$(CXX) $(CPPFLAGS) -c -o rand_serial.o rand_serial.cc 

admixmap:	$(infiles)
		$(CXX) $(CPPFLAGS) -o ../test/admixmap $(infiles) $(LIBS)

check:		all
		perl ../test/batchtest.pl

include: .depend

#Compile  for Linux and repackage distribution with new exec
linuxdist:      all
		cp ../test/admixmap $(DISTDIR)
		rm -f ../test/$(DISTDIR).tar.gz
		tar -czf admixmap.tar.gz ../test/$(DISTDIR)

#Compile for Windows and prepare to repackage
#Zipping cannot be done in Windows from command line unless a special zip app installed (eg Winzip, AlZip etc)
windist:	all
		cp ../test/admixmap.exe $(WINDISTDIR)
		rm -f ../test/$(WINDISTDIR).zip
## ALZip
		alzip -a ../test/$(WINDISTDIR) ../test/admixmap.zip
## WinZip
#		??
## pkzip
# 		pkzip admixmap.zip ../test/$(WINDISTDIR)


