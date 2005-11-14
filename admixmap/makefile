# makefile for compiling ADMIXMAP 

CC			= gcc
CXX 			= g++
DISTDIR			= ../test/admixmap-linux
WINDISTDIR		= ../test/admixmap-win
LIBS			= -lm -static -lgsl -lgslcblas
# processor architecture flags
I686FLAGS    = -march=i686 
P4FLAGS      = -march=pentium4 -mfpmath=sse
X86_64FLAGS  = -march=athlon64
# flags for release version, debug version or profiling version
DFLAGS 		= -O3		# release
# DFLAGS 	= -g    	# debug
# DFLAGS 	= -pg -O3 	# profiling

CXXFLAGS	= 
CPPFLAGS	= -W -Wall $(I686FLAGS) $(DFLAGS)
# CPPFLAGS	= -W -Wall $(P4FLAGS) $(DFLAGS)

objects		= rand_serial.o admixmap.o Latent.o Regression.o AlleleFreqs.o AdmixOptions.o \
AdaptiveRejection.o StepSizeTuner.o Gaussian.o Genome.o Chromosome.o CompositeLocus.o IndividualCollection.o \
Individual.o GaussianProposalMH.o HMM.o IndAdmixOutputter.o functions.o chib.o StratificationTest.o ScoreTests.o \
MisSpecAlleleFreqTest.o DispersionTest.o HWTest.o LogWriter.o InputData.o StringSplitter.o StringConvertor.o \
DirichletParamSampler.o HamiltonianMonteCarlo.o DataMatrix.o MuSampler.o DispersionSampler.o

all:		admixmap .depend

again: 		clean all

clean:		          
		rm -f *.o
		rm -f .depend

.depend: 	*.cc
		$(CXX) $(CPPFLAGS) -MM *.cc >> .depend # -MM option outputs rule describing dependencies 

admixmap:	$(objects)
		$(CXX) $(CPPFLAGS) -o ../test/admixmap $(objects) $(LIBS)

check:		all
		perl ../test/batchtest.pl

include: 	.depend 

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


