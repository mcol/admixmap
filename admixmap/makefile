# makefile for compiling ADMIXMAP
# edit to suit your system by commenting/uncommenting lines (with #)

##default is to make serial version with gcc compiler. To make parallel version, edit lines below or use make VERSION=parallel
##remember to use make [VERSION=...] new when switching between versions
##for debug (or profiling), use make [VERSION=...] DFLAGS=-g (or-pg) new

VERSION = serial
#VERSION = parallel

##paths
GSL_PATH =#only needed if in nonstandard location 
MPICH_PATH = /usr/local/mpich/gcc##parallel only
#MPICH_PATH = /usr/local/mpich2/gcc##parallel only, mpich2
SPRNG_PATH = /opt/packages/sprng-2.0##parallel only
MPICH_LIBDIR = lib64#this sometimes has different names

## ** Choose your preferred compiler
GNU_COMPILER  = g++ # GNU compiler, all machines
INTEL_COMPILER = icpc# Intel compiler for Hamilton
PORTLAND_COMPILER = pgCC# Portland Group compiler, for Walton
MPICC = $(MPICH_PATH)/bin/mpiCC# MPI wrapper for gnu, MPI 1
MPICXX = $(MPICH_PATH)/bin/mpicxx# MPI wrapper for gnu, MPI 2

## **Destination details
DESTDIR = ../test# where to put compiled exec
SERIAL_EXEC = admixmap# name of compiled executable
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

## ** 
CPPFLAGS = $(GNUFLAGS)
#CPPFLAGS = $(ICPCFLAGS)
#CPPFLAGS = $(MPICCFLAGS)
#CPPFLAGS = $(PGCCFLAGS)

## Library flags
# math, GSL, GSL's CBLAS (linear algebra)
SERIAL_LIBS = -lm -lgsl -lgslcblas  #-static # serial version
##NB: do not use -static on cluster

PARALLEL_LIBS = -lm -lgsl -lgslcblas -lmpich -lsprng -lgmp -lmpe# parallel version
#mpich, sprng(parallel RNG), GNU message passing

##library flags for sprng in parallel version
PARALLEL_LFLAGS = -L$(SPRNG_PATH)/lib -L$(MPICH_PATH)/$(MPICH_LIBDIR)  # -L/usr/local/lib
PARALLEL_INCLUDES = -I$(SPRNG_PATH)/include  -I$(MPICH_PATH)/include# -I/usr/local/include

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
EXEC = $(SERIAL_EXEC)
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

new: 		clean all

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


