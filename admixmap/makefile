# makefile for compiling ADMIXMAP
# edit to suit your system by commenting/uncommenting lines (with #)

#VERSION = serial
VERSION = parallel

## ** Choose your preferred compiler
GNU_COMPILER  = g++ # GNU compiler, all machines
INTEL_COMPILER = icpc# Intel compiler for Hamilton
PORTLAND_COMPILER = pgCC# Portland Group compiler, for Walton
MPICC = /usr/local/mpich/gcc/bin/mpiCC# MPI wrapper for gnu, MPI 1
#CXX = /usr/local/mpich2/gcc/bin/mpicxx# MPI wrapper for gnu, MPI 2

## **Destination details
DESTDIR = ../test# where to put compiled exec
#EXEC = adm-para# name of compiled executable
DISTDIR			= $(DESTDIR)/admixmap-linux#directory for Linux distribution
WINDISTDIR		= $(DESTDIR)/admixmap-win#  directory for Windows distribution

## processor architecture flags
## ** change this to your processor type
## g++
#PROC = -march=i686  # Intel 686
#PROC = -march=pentium4 -mfpmath=sse # Pentium4
PROC = -march=k8  # Walton, g++
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
MPICCFLAGS = -w0 $(DFLAGS)#            mpiCC, ?redundant, use GNUFLAGS
PGCCFLAGS = $(PROC) $(DFLAGS)#         pgCC

## Library flags
# math, GSL, GSL's CBLAS (linear algebra)
SERIAL_LIBS = -lm -lgsl -lgslcblas  #-static # serial version
PARALLEL_LIBS = -lm -lgsl -lgslcblas -lmpich -lsprng -lgmp -lmpe# parallel version
#mpich, sprng(parallel RNG), GNU message passing
##NB: do not use -static on cluster

##library flags for sprng in parallel version
PARALLEL_LFLAGS = -L/opt/packages/sprng-2.0/lib -L/usr/local/mpich/gcc/lib64  # -L/usr/local/lib
PARALLEL_INCLUDES = -I/opt/packages/sprng-2.0/include  -I/usr/local/mpich/gcc/include# -I/usr/local/include

## ** 
CPPFLAGS = $(GNUFLAGS)
#CPPFLAGS = $(ICPCFLAGS)
#CPPFLAGS = $(MPICCFLAGS)
#CPPFLAGS = $(PGCCFLAGS)

ifeq ($(VERSION),parallel)
CXX = $(MPICC)
LIBS = $(PARALLEL_LIBS)
LFLAGS = $(PARALLEL_LFLAGS)
INCLUDES = $(PARALLEL_INCLUDES)
EXEC = adm-para
else
CXX = $(GNU_COMPILER)
LIBS = $(SERIAL_LIBS)
LFLAGS = 
INCLUDES = 
EXEC = admixmap
endif

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


