###############################################################################################
# makefile for compiling ADMIXMAP
#
# The default is to make serial version with gcc compiler and optimised for i686 machines. 
# For other versions, edit the following variables or specify on the command line. 
# Use make[...] new to recompile from scratch.
#
# Serial/Parallel version: VERSION = [serial/parallel]
# Processor type: ARCH = [i686/P4/AMD64/itanium2]
# Compiler: CC = [g++/icpc/pgCC/pathCC]
# MPI wrapper for compiler PCC = [mpiCC/mpicxx]
# debug/profiling/optimisation : DFLAGS = [-g/-pg/-O3] (NB: do not optimise with -g)
# paths to gsl, mpi, sprng: GSL_LIBPATH, GSL_INCLUDEPATH, MPI_PATH, MPI_LIBPATH, SPRNG_PATH
# name of compiled executable: EXEC
#
# [the following does not apply to Windows]
# To use shared libraries, leave LIBTYPE blank and add the paths to gsl 
# (and possibly mpi and sprng in parallel version) shared libs (.so) to your dynamic library path (LD_LIBRARY_PATH in bash)
# To use static libraries, specify LIBTYPE=-static.
#
###############################################################################################

# GSL - only needed if in nonstandard location 
GSL_LIBPATH =#standard is /usr/local/lib
GSL_INCLUDEPATH=#standard is /usr/local/include

## for parallel version only
# path to mpi top level directory
MPI_PATH = /usr/local/mpich/gcc# mpich with gcc 3
#MPI_PATH = /usr/local/mpich/gcc4# mpich with gcc 4
#MPI_PATH = /usr/local/mpich/path## MPI wrapper for Pathscale compiler
#MPI_PATH = /usr/local/mpich2/gcc## mpich2 on Walton
MPI_LIBPATH = $(MPIPATH)/lib64# this sometimes has different names

# path to sprng top level directory
SPRNG_PATH = /opt/packages/sprng-2.0

#leave this one empty for shared libraries (eg libgsl.so) or as -static for static libs (eg libgsl.a)
#NB do not use -static with parallel version
LIBTYPE =-static

## compiler definitions
GNU_COMPILER  = g++ # GNU compiler, all machines
INTEL_COMPILER = icpc# Intel compiler for Hamilton
PATHSCALE_COMPILER = pathCC# Pathscale compiler,( Walton )
PORTLAND_COMPILER = pgCC# Portland Group compiler ( Walton )
MPICC = $(MPI_PATH)/bin/mpiCC# MPI 1 wrapper
MPICXX = $(MPI_PATH)/bin/mpicxx# MPI 2 wrapper

## DEFAULTS
VERSION = serial
#VERSION = parallel

## processor type
ARCH = i686# Intel 686
#ARCH = P4#Pentium4
#ARCH = AMD64# Walton
#ARCH= itanium2# Hamilton (default)

# compiler
CC = $(GNU_COMPILER)#serial
PCC = $(MPICC)#parallel wrapper for serial

## **Destination details
DESTDIR = ../test# where to put compiled exec
SERIAL_EXEC = admixmap# name of compiled executable

## ** flags for release version, debug version or profiling version
DFLAGS 		= -O3# release, optimized
#DFLAGS 	= -g# debug
# DFLAGS 	= -pg -O3# profiling

#warning flags for compiler
WFLAGS =

# processor architecture flags
PROC =

## DO NOT EDIT BELOW HERE, except to add processor/compiler types
ifeq ($(CC),$(GNU_COMPILER))
WFLAGS = -W -Wall
ifeq ($(ARCH),i686)
PROC = -march=i686  # Intel 686
else
ifeq ($(ARCH),P4)
PROC = -march=pentium4 -mfpmath=sse # Pentium4
else
ifeq ($(ARCH),AMD64)
PROC = -march=k8  
endif
endif
endif
endif

ifeq ($(CC),$(INTEL_COMPILER))
WFLAGS = -w0
ifeq ($(ARCH),itanium2)
PROC = -mcpu=itanium2
endif
endif

ifeq ($(CC),$(PATHSCALE_COMPILER))
ifeq ($(ARCH),AMD64)
PROC = -mcpu=auto
endif
endif

ifeq ($(CC),$(PORTLAND_COMPILER))
ifeq ($(ARCH),AMD64)
PROC = -tp amd64
#DFLAGS = $(DFLAGS) -ipa
endif
endif

####  DO NOT EDIT BELOW HERE ####

## Library flags
# math, GSL, GSL's CBLAS (linear algebra)
SERIAL_LIBS = -lm -lgsl -lgslcblas 

PARALLEL_LIBS = $(SERIAL_LIBS) -lmpich -lsprng -lgmp -lmpe
#mpich, sprng(parallel RNG), GNU message passing

ifeq ($(GSL_LIBPATH),)
GSL_LFLAGS =
else
GSL_LFLAGS = -L$(GSL_LIBPATH)
endif

ifeq ($(GSL_INCLUDEPATH),)
GSL_INCLUDES =
else
GSL_INCLUDES = -I$(GSL_INCLUDEPATH)
endif

##library flags for sprng and mpich in parallel version
PARALLEL_LFLAGS = -L$(SPRNG_PATH)/lib -L$(MPI_LIBPATH)
PARALLEL_INCLUDES = -I$(SPRNG_PATH)/include  -I$(MPI_PATH)/include

## parallel version
ifeq ($(VERSION),parallel)
CXX = $(PCC)
LIBS = $(PARALLEL_LIBS)
LFLAGS = $(GSL_LFLAGS) $(PARALLEL_LFLAGS)
INCLUDES = $(GSL_INCLUDES) $(PARALLEL_INCLUDES)
EXEC = adm-para
##serial version
else
CXX = $(CC)
LIBS = $(SERIAL_LIBS) $(LIBTYPE)
LFLAGS = $(GSL_LFLAGS)
INCLUDES = $(GSL_INCLUDES)
EXEC = $(SERIAL_EXEC)
endif

## **compiler/linker flags
CPPFLAGS = $(WFLAGS) $(PROC) $(DFLAGS)

objects		= rand.o admixmap.o Latent.o Regression.o AlleleFreqs.o AdmixOptions.o \
AdaptiveRejection.o StepSizeTuner.o Gaussian.o Genome.o Chromosome.o CompositeLocus.o IndividualCollection.o \
Individual.o GaussianProposalMH.o HMM.o IndAdmixOutputter.o functions.o chib.o StratificationTest.o ScoreTests.o \
ScoreTestBase.o MisSpecAlleleFreqTest.o DispersionTest.o HWTest.o LogWriter.o InputData.o StringSplitter.o StringConvertor.o \
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

