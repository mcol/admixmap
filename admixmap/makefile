# makefile for compiling ADMIXMAP

CXX  = g++ # GNU compiler, all machines
#CXX = icpc# Intel compiler for Hamilton
#CXX = pgCC# Portland Group compiler, for Walton

DISTDIR			= ../test/admixmap-linux
WINDISTDIR		= ../test/admixmap-win
LIBS			= -lm -static -lgsl -lgslcblas
EXEC = ../test/admixmap #where to put compiled exec

#  g++ processor architecture flags
I686FLAGS    = -march=i686 
P4FLAGS      = -march=pentium4 -mfpmath=sse
X86_64FLAGS  = -march=athlon64
K8FLAGS      = -march=k8#Walton
ITANIUMFLAGS = -mcpu=itanium2#Hamilton 

ifeq ($(CXX), g++)
LFLAGS = #path to extra libs
INCLUDES = #path to extra header files
PROC = $(I686FLAGS)
WFLAGS = -W -Wall #Warning flags
endif 
ifeq ($(CXX),pgCC)#Walton Flags
PROJDIR = /ichec/work/gsl#path to gsl
LFLAGS = -L$(PROJDIR)/lib
INCLUDES = -I$(PROJDIR)/include
PROC = -tp amd64
WFLAGS = 
endif
ifeq ($(CXX),icpc)#Hamilton flags
PROJDIR = /ichec/home/users/doducd/share#path to gsl
LFLAGS = -L$(PROJDIR)/lib
INCLUDES = -I$(PROJDIR)/include
WFLAGS = -w0#w0 for errors but no warnings, w1 for warnings too, w2 for comments
PROC = #-mcpu=itanium2 #default option
endif

# flags for release version, debug version or profiling version
DFLAGS 		= -O3# release
#DFLAGS 	= -g# debug
# DFLAGS 	= -pg -O3# profiling

CXXFLAGS	= 
CPPFLAGS	= $(WFLAGS) $(PROC) $(DFLAGS)

objects		= rand_serial.o admixmap.o Latent.o Regression.o AlleleFreqs.o AdmixOptions.o \
AdaptiveRejection.o StepSizeTuner.o Gaussian.o Genome.o Chromosome.o CompositeLocus.o IndividualCollection.o \
Individual.o GaussianProposalMH.o HMM.o IndAdmixOutputter.o functions.o chib.o StratificationTest.o ScoreTests.o \
MisSpecAlleleFreqTest.o DispersionTest.o HWTest.o LogWriter.o InputData.o StringSplitter.o StringConvertor.o \
DirichletParamSampler.o HamiltonianMonteCarlo.o DataMatrix.o MuSampler.o DispersionSampler.o

all:  	$(objects)
	$(CXX) $(CPPFLAGS) -o $(EXEC) $(objects) $(LFLAGS) $(LIBS)

%.o: %.cc#rule for compilation
	$(CXX) $(CPPFLAGS) $(INCLUDES) -c $< -o $@

again: 		clean all

clean:		          
	rm -f *.o

realclean: 
	rm -f *.o 
	rm -f $(EXEC)

check:		all
		perl ../test/batchtest.pl

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


