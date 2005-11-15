# makefile for compiling ADMIXMAP 

CC			= gcc
CXX 			= g++
DISTDIR			= ../test/admixmap-linux
WINDISTDIR		= ../test/admixmap-win
LIBS			= -lm -static -lgsl -lgslcblas
EXEC = ../test/admixmap #where to put compiled exec

LFLAGS = #path to extra libs
INCLUDES = #path to extra header files

#Walton Flags
#PROJDIR = /ichec/work/ndlif004
#LFLAGS = -L$(PROJDIR)/lib
#INCLUDES = -I$(PROJDIR)/include

# processor architecture flags
I686FLAGS    = -march=i686 
P4FLAGS      = -march=pentium4 -mfpmath=sse
X86_64FLAGS  = -march=athlon64
K8FLAGS      = -march=k8 

# flags for release version, debug version or profiling version
DFLAGS 		= -O3		# release
# DFLAGS 	= -g    	# debug
# DFLAGS 	= -pg -O3 	# profiling

PROC = $(I686FLAGS)
#PROC = $(K8FLAGS)##for Walton
CXXFLAGS	= 
CPPFLAGS	= -W -Wall $(PROC) $(DFLAGS)

objects		= rand_serial.o admixmap.o Latent.o Regression.o AlleleFreqs.o AdmixOptions.o \
AdaptiveRejection.o StepSizeTuner.o Gaussian.o Genome.o Chromosome.o CompositeLocus.o IndividualCollection.o \
Individual.o GaussianProposalMH.o HMM.o IndAdmixOutputter.o functions.o chib.o StratificationTest.o ScoreTests.o \
MisSpecAlleleFreqTest.o DispersionTest.o HWTest.o LogWriter.o InputData.o StringSplitter.o StringConvertor.o \
DirichletParamSampler.o HamiltonianMonteCarlo.o DataMatrix.o MuSampler.o DispersionSampler.o

all:  	$(objects)
	$(CXX) $(CPPFLAGS) -o $(EXEC) $(objects) $(LFLAGS) $(LIBS)

%.o: %.cc
	$(CXX) $(CPPFLAGS) $(INCLUDES) -c $< -o $@

again: 		clean all

clean:		          
	rm -f *.o

realclean: 
	rm -f *.o 
	rm -f $(EXEC)

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


