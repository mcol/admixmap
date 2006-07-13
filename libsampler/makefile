###############################################################################################
# makefile for compiling libsampler, the sampling library
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
# paths to gsl, mpi, sprng headers: GSL_INCLUDEPATH, MPI_PATH, SPRNG_PATH
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

## compiler definitions
GNU_COMPILER  = g++ # GNU compiler, all machines
INTEL_COMPILER = icpc# Intel compiler for Hamilton
PATHSCALE_COMPILER = pathCC# Pathscale compiler,( Walton )
PORTLAND_COMPILER = pgCC# Portland Group compiler ( Walton )
MPICC = $(MPI_PATH)/bin/mpiCC# MPI 1 wrapper
MPICXX = $(MPI_PATH)/bin/mpicxx# MPI 2 wrapper

## DEFAULTS
GCC_VERSION =3
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
DFLAGS=-O3 -ftree-vectorize
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
WFLAGS = -fullwarn -LNO:simd_verbose=ON
#WFLAGS = -Wall 
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

####  DO NOT EDIT BELOW HERE ####

ifeq ($(GSL_INCLUDEPATH),)
GSL_INCLUDES =
else
GSL_INCLUDES = -I$(GSL_INCLUDEPATH)
endif

##include flags for sprng and mpich in parallel version
PARALLEL_INCLUDES = -I$(SPRNG_PATH)/include  -I$(MPI_PATH)/include

## parallel version
ifeq ($(VERSION),parallel)
CXX = $(PCC)
INCLUDES = $(GSL_INCLUDES) $(PARALLEL_INCLUDES)
##serial version
else
CXX = $(CC)
INCLUDES = $(GSL_INCLUDES)
endif

## **compiler/linker flags
CPPFLAGS = $(WFLAGS) $(PROC) $(DFLAGS)
PICFLAG =# -fPIC

## DO NOT EDIT BELOW HERE ##

objects		=  AdaptiveRejection.o StepSizeTuner.o Gaussian.o GaussianProposalMH.o \
DirichletParamSampler.o HamiltonianMonteCarlo.o MuSampler.o DispersionSampler.o rand.o functions.o

all:	static #shared

static: libsampler

libsampler: $(objects) 
	ar crs libsampler.a $(objects)
	#ranlib libsampler.a

shared: $(objects)
	$(CXX) $(CPPFLAGS) -shared -Wl,-soname,libsampler.so\
-o libsampler.so $(objects) -lgsl -lgslcblas

%.o: %.cc#rule for compilation
	$(CXX) $(CPPFLAGS) $(INCLUDES) $(PICFLAG) -c $< -o $@

clean:		          
	rm -f *.o

veryclean: clean
	rm -f *.a *.so

dist:
	rm libsampler.tar.gz
	mkdir libsampler
	cp *.cc *.h *.txt README libsampler
	cp makefile libsampler
	tar -czf libsampler.tar.gz libsampler
	rm -R libsampler




