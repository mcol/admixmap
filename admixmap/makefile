###############################################################################################
# makefile for compiling ADMIXMAP
#
# The default is to make serial version with gcc compiler and optimised for i686 machines. 
# For other versions, edit the following variables or specify on the command line.
# Use make bayeslib to compile the sampling library if not already done.
# Use make[...] new to recompile from scratch or make[...] allnew to recompile the sampling library as well.
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
GSL_LIBPATH =/usr/local/lib
GSL_INCLUDEPATH=/usr/local/include

#path to sampling library
BAYESLIB_PATH =../bayeslib

## for parallel version only
# path to mpi top level directory
#MPI_PATH = /usr/local/mpich/gcc# mpich with gcc 3
#MPI_PATH = /usr/local/mpich/gcc4# mpich with gcc 4
MPI_PATH = /usr/local/mpich/path## MPI wrapper for Pathscale compiler
#MPI_PATH = /usr/local/mpich2/gcc## mpich2 on Walton
MPI_LIBPATH = $(MPI_PATH)/lib64# this sometimes has different names

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
GCC_VERSION =4
VERSION = serial
#VERSION = parallel

## processor type
#ARCH = i686# Intel 686
#ARCH = P4#Pentium4
#ARCH = AMD64# Walton
ARCH= itanium2# Hamilton (default)

# compiler
CC = $(INTEL_COMPILER)#serial
PCC = $(MPICC)#parallel wrapper for serial

## **Destination details
DESTDIR = ../test# where to put compiled exec
EXEC = admixmap# name of compiled executable

## ** flags for release version, debug version or profiling version
DFLAGS 		= -O3# release, optimized
#DFLAGS 	= -g# debug
# DFLAGS 	= -pg -O3# profiling


#warning flags for compiler
WFLAGS =

# processor architecture flags
PROC =

#GCC_VERSION = $(shell gcc -dumpversion)
## DO NOT EDIT BELOW HERE, except to add processor/compiler types
ifeq ($(CC),$(GNU_COMPILER))
WFLAGS = -W -Wall
ifeq ($(GCC_VERSION),4)
ifeq ($(DFLAGS),-O3)
override DFLAGS+= -ftree-vectorize
endif	
endif
ifeq ($(ARCH),i686)
PROC = -march=i686# Intel 686
else
ifeq ($(ARCH),P4)
PROC = -march=pentium4 -msse2 -mfpmath=sse# Pentium4
else
ifeq ($(ARCH),AMD64)
PROC = -march=k8 -msse2 -mfpmath=sse 
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
#WFLAGS = -fullwarn -LNO:simd_verbose=ON#for info on vectorization
WFLAGS = -Wall 
ifeq ($(ARCH),AMD64)
PROC = -mcpu=auto
DFLAGS = -O3 -OPT:Ofast -fno-math-errno -ffast-math #equivalent to -Ofast but without ipa
endif
endif

ifeq ($(CC),$(PORTLAND_COMPILER))
ifeq ($(ARCH),AMD64)
PROC = -tp amd64
endif
endif

####  DO NOT EDIT BELOW HERE ##

## Library flags
# math, GSL, GSL's CBLAS (linear algebra)
NORMAL_LIBS = $(BAYESLIB_PATH)/bayeslib.a -lm -lgsl -lgslcblas 

PARALLEL_LIBS = $(NORMAL_LIBS) -lmpich -lsprng -lgmp -lmpe
#mpich, sprng(parallel RNG), GNU message passing

LFLAGS =
INCLUDES =

ifneq ($(GSL_LIBPATH),)
LFLAGS += -L$(GSL_LIBPATH)
endif

ifneq ($(GSL_INCLUDEPATH),)
INCLUDES += -I$(GSL_INCLUDEPATH)
endif

##library flags for sprng and mpich in parallel version
PARALLEL_LFLAGS = -L$(SPRNG_PATH)/lib -L$(MPI_LIBPATH)
PARALLEL_INCLUDES = $(INCLUDES) -I$(SPRNG_PATH)/include  -I$(MPI_PATH)/include
PARALLEL_EXEC = adm-para
SERIAL_INCLUDES=$(INCLUDES)

LIBS = $(LIBTYPE) $(NORMAL_LIBS)
CXX=$(CC)

## **compiler/linker flags
CPPFLAGS = $(WFLAGS) $(PROC) $(DFLAGS)

export CC PCC CPPFLAGS SERIAL_INCLUDES PARALLEL_INCLUDES

objects	= admixmap.o Latent.o AlleleFreqs.o AdmixOptions.o Genome.o Chromosome.o CompositeLocus.o\
IndividualCollection.o Individual.o HMM.o IndAdmixOutputter.o chib.o StratificationTest.o ScoreTestBase.o ScoreTests.o \
AdmixtureAssocTest.o AffectedsOnlyTest.o ResidualLDTest.o MisSpecAlleleFreqTest.o DispersionTest.o HWTest.o InputData.o AlleleFreqSampler.o Comms.o

all:	checkbayeslib message $(objects)
	$(CXX) $(CPPFLAGS) -o $(DESTDIR)/$(EXEC) $(objects) $(LFLAGS) $(LIBS)
	@echo **ADMIXMAP has been compiled as $(DESTDIR)/$(EXEC) \**

serial: all

parallel:
	@$(MAKE) INCLUDES="$(PARALLEL_INCLUDES)" LFLAGS="$(PARALLEL_LFLAGS)" LIBS="$(PARALLEL_LIBS)" CXX=$(PCC) EXEC=$(PARALLEL_EXEC) BAYESLIB_RULES=parallel

message:	
	@echo **Compiling ADMIXMAP**

checkbayeslib:
	test -e $(BAYESLIB_PATH)/bayeslib.a || $(MAKE) -C $(BAYESLIB_PATH) $(BAYESLIB_RULES) -e

bayeslib:
	$(MAKE) -C $(BAYESLIB_PATH) -e
samplers:
	@echo **compiling sampler library**
	$(MAKE) -C $(BAYESLIB_PATH)/samplers -e

utils: 
	@echo **compiling utility library**
	$(MAKE) -C $(BAYESLIB_PATH)/utils -e


%.o: %.cc#rule for compilation
	$(CXX) $(CPPFLAGS) -I$(BAYESLIB_PATH) $(INCLUDES) -c $< -o $@

#recompile from scratch
new: 		clean all

#recompile everything from scratch
allnew:	libnew new

#delete all object files
clean:		          
	rm -f *.o

#delete sampler library
libclean: 
	$(MAKE) -C $(BAYESLIB_PATH) clean

#recompile sampler library
libnew: 
	@echo **recompiling bayeslib**
	$(MAKE) -C $(BAYESLIB_PATH) -e allnew bayeslib

allclean: libclean clean

#delete all object files, libraries and exec
realclean: clean libclean 
	rm -f $(DESTDIR)/$(EXEC)

check:		all
		perl $(DESTDIR)/batchtest.pl




