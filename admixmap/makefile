# makefile for compiling ADMIXMAP
# edit to suit your system by commenting/uncommenting lines (with #)

## ** Choose your preferred compiler
CXX  = g++ # GNU compiler, all machines
#CXX = icpc# Intel compiler for Hamilton
#CXX = pgCC# Portland Group compiler, for Walton
#CXX = mpiCC#

## **Destination details
DESTDIR = ../test# where to put compiled exec
EXEC = admixmap# name of compiled executable
DISTDIR			= $(DESTDIR)/admixmap-linux#directory for Linux distribution
WINDISTDIR		= $(DESTDIR)/admixmap-win#  directory for Windows distribution

## processor architecture flags
## ** change this to your processor type
## g++
#PROC = -march=i686  # Intel 686
PROC = -march=pentium4 -mfpmath=sse # Pentium4
#PROC = -march=k8  # Walton, g++
## icpc
#PROC = -mcpu=itanium2  #Hamilton (default)
## pgCC
#PROC = -tp amd64  # Walton

## ** flags for release version, debug version or profiling version
DFLAGS 		= -O3# release, optimized
#DFLAGS 	= -g# debug
# DFLAGS 	= -pg -O3# profiling

## **compiler/linker flags
GNUFLAGS = -W -Wall $(PROC) $(DFLAGS)# g++ compiler
ICPCFLAGS = -w0 $(DFLAGS)#             icpc
MPICCFLAGS = -w0 $(DFLAGS)#            mpiCC
PGCCFLAGS = $(PROC) $(DFLAGS)#         pgCC

## Library flags
LIBS = -lm -static -lgsl -lgslcblas  # serial version
#LIBS = -lm -lgsl -lgslcblas -lmpich # parallel version
LFLAGS =                             # -L/usr/local/lib
INCLUDES =                           # -I/usr/local/include

## ** 
CPPFLAGS = $(GNUFLAGS)
#CPPFLAGS = $(ICPCFLAGS)
#CPPFLAGS = $(MPICCFLAGS)
#CPPFLAGS = $(PGCCFLAGS)


## DO NOT EDIT BELOW HERE ##

objects		= rand.o admixmap.o Latent.o Regression.o AlleleFreqs.o AdmixOptions.o \
AdaptiveRejection.o StepSizeTuner.o Gaussian.o Genome.o Chromosome.o CompositeLocus.o IndividualCollection.o \
Individual.o GaussianProposalMH.o HMM.o IndAdmixOutputter.o functions.o chib.o StratificationTest.o ScoreTests.o \
MisSpecAlleleFreqTest.o DispersionTest.o HWTest.o LogWriter.o InputData.o StringSplitter.o StringConvertor.o \
DirichletParamSampler.o HamiltonianMonteCarlo.o DataMatrix.o MuSampler.o DispersionSampler.o AlleleFreqSampler.o

all:  	$(objects)
	$(CXX) $(CPPFLAGS) -o $(DESTDIR)/$(EXEC) $(objects) $(LFLAGS) $(LIBS)

%.o: %.cc#rule for compilation
	$(CXX) $(CPPFLAGS) $(INCLUDES) -c $< -o $@

again: 		clean all

parallel:#TODO 

clean:		          
	rm -f *.o

realclean: 
	rm -f *.o 
	rm -f $(DESTDIR)/$(EXEC)

check:		all
		perl $(DESTDIR)/batchtest.pl

#Compile  for Linux and repackage distribution with new exec
linuxdist:      all
		cp $(DESTDIR)/$(EXEC) $(DISTDIR)
		rm -f $(DESTDIR)/$(EXEC).tar.gz
		tar -czf $(EXEC).tar.gz $(DISTDIR)

#Compile for Windows and prepare to repackage
#Zipping cannot be done in Windows from command line unless a special zip app installed (eg Winzip, AlZip etc)
windist:	all
		cp $(DESTDIR)/$(EXEC).exe $(WINDISTDIR)
		rm -f $(DESTDIR)/$(EXEC).zip
## ALZip
		alzip -a $(WINDISTDIR) $(DESTDIR)/$(EXEC).zip
## WinZip
#		??
## pkzip
# 		pkzip $(EXEC).zip $(WINDISTDIR)


