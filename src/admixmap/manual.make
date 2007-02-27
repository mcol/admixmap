###############################################################################################
# makefile for compiling ADMIXMAP and HAPMIXMAP
#
# The default is to make serial version with gcc compiler and optimised for i686 machines. 
# For other versions, edit the variables listed below or specify on the command line.
# Use 'make admixmap or make hapmixmap for separate compilation. Use 'make' or 'make all' for both.
# Use 'make bayeslib' to compile the sampling library if not already done.
# Use 'make clean' to remove all object files, 'make allclean' also removes bayeslib.
# Use 'make new' to recompile from scratch or 'make allnew' to also recompile the sampling library.
#
# Processor type: ARCH = [i686/P4/AMD64/itanium2]
# Compiler: CC = [g++/icpc/pgCC/pathCC]
# (gcc compiler only) GCC_VERSION =3/4 (NB no space after =)
# MPI wrapper for compiler PCC = [mpiCC/mpicxx]
# debug/profiling/optimisation : DFLAGS = [-g/-pg/-O3] (NB: do not optimise with -g)
# paths to gsl, mpi, sprng: GSL_LIBPATH, GSL_INCLUDEPATH, MPI_PATH, MPI_LIBPATH, SPRNG_PATH
# name of compiled executable: ADMEXEC, HAPEXEC
#
# [the following does not apply to Windows]
# To use shared libraries, leave LIBTYPE blank and add the paths to gsl 
# (and possibly mpi and sprng in parallel version) shared libs (.so) to your dynamic library path (LD_LIBRARY_PATH in bash)
# To use static libraries, specify LIBTYPE=-static.
#
###############################################################################################

# GSL - only needed if in nonstandard location 
GSL_LIBPATH =#/usr/local/lib
GSL_INCLUDEPATH=#/usr/local/include

#path to sampling library
BAYESLIB_PATH =../bayeslib
BAYESLIB_NAME =bayeslib.a

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
ARCH = AMD64# Walton
#ARCH= itanium2# Hamilton (default)

# compiler
CC = $(PATHSCALE_COMPILER)#serial
#CC = $(GNU_COMPILER)
PCC = $(MPICC)#parallel wrapper for serial

## **Destination details
DESTDIR = ../test# where to put compiled exec
ADMEXEC = admixmap# name of compiled executable
HAPEXEC = hapmixmap

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
DFLAGS = -O3 -OPT:Ofast -fno-math-errno -ffast-math -fexceptions#equivalent to -Ofast but without ipa
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
NORMAL_LIBS = -lm -lgsl -lgslcblas 

PARALLEL_LIBS =  $(BAYESLIB_PATH)/parabayeslib.a $(NORMAL_LIBS) -lmpich -lsprng -lgmp -lmpe
#mpich, sprng(parallel RNG), GNU message passing

LFLAGS =#-L$(SPRNG_PATH)/lib
INCLUDES =

ifneq ($(GSL_LIBPATH),)
LFLAGS += -L$(GSL_LIBPATH)
endif

ifneq ($(GSL_INCLUDEPATH),)
INCLUDES += -I$(GSL_INCLUDEPATH)
endif

##library flags for sprng and mpich in parallel version
PARALLEL_LFLAGS = -L$(SPRNG_PATH)/lib -L$(MPI_LIBPATH)
PARALLEL_INCLUDES = $(INCLUDES) -DPARALLEL -I$(SPRNG_PATH)/include  -I$(MPI_PATH)/include
PARALLEL_ADMEXEC = admixmap-para
PARALLEL_HAPEXEC = hapmixmap-para
SERIAL_INCLUDES=$(INCLUDES) 
PARALLEL_DEFINES=

LIBS = $(LIBTYPE) $(BAYESLIB_PATH)/bayeslib.a $(NORMAL_LIBS)
CXX=$(CC)

## **compiler/linker flags
CPPFLAGS = $(WFLAGS) $(PROC) $(DFLAGS)

export CC PCC CPPFLAGS SERIAL_INCLUDES PARALLEL_INCLUDES

common_objects	= misc.o Model.o AlleleFreqs.o Options.o AdmixOptions.o InputData.o Genome.o\
 Chromosome.o CompositeLocus.o Haplotype.o IndividualCollection.o Individual.o HMM.o AdmixtureAssocTest.o\
 AffectedsOnlyTest.o AncestryAssocTest.o ScoreTestBase.o ScoreTests.o \
 ResidualLDTest.o HWTest.o AlleleFreqSampler.o Comms.o GenotypeProbOutputter.o Annealer.o

admixmap_objects = admixmap.o AdmixMapModel.o PopAdmix.o DispersionFreqs.o AdmixedIndividual.o AdmixIndividualCollection.o\
 chib.o IndAdmixOutputter.o StratificationTest.o MisSpecAlleleFreqTest.o DispersionTest.o

hapmixmap_objects = hapmixmap.o HapMixModel.o HapMixOptions.o PopHapMix.o\
HapMixFreqs.o HapMixIndividual.o HapMixIndividualCollection.o MantelHaenszelTest.o 

#both serial
all:	admixmap hapmixmap

#ADMIXMAP serial
admixmap:	checkbayeslib admixmap_message define_admixmap $(admixmap_objects) $(common_objects)
	$(CXX) $(CPPFLAGS) -o $(DESTDIR)/$(ADMEXEC) $(admixmap_objects) $(common_objects) $(LFLAGS) $(LIBS)
	@echo **ADMIXMAP has been compiled as $(DESTDIR)/$(ADMEXEC) \**
#HAPMIXMAP serial
hapmixmap:	checkbayeslib hapmixmap_message define_hapmixmap $(hapmixmap_objects) $(common_objects)
	$(CXX) $(CPPFLAGS) -o $(DESTDIR)/$(HAPEXEC) $(hapmixmap_objects) $(common_objects) $(LFLAGS) $(LIBS)
	@echo **HAPMIXMAP has been compiled as $(DESTDIR)/$(HAPEXEC) \**

serial: all

#NOTE defaults to parallel version of hapmixmap
parallel: hapmixmap-para

#HAPMIXMAP parallel
hapmixmap-para:
	@$(MAKE) -fmanual.make hapmixmap INCLUDES="$(PARALLEL_INCLUDES)" LFLAGS="$(PARALLEL_LFLAGS)" LIBS="$(PARALLEL_LIBS)" CXX=$(PCC) HAPEXEC=$(PARALLEL_HAPEXEC) BAYESLIB_RULES=parallel BAYESLIB_NAME=parabayeslib.a
#ADMIXMAP parallel
admixmap-para:
	@$(MAKE) -fmanual.make admmixmap INCLUDES="$(PARALLEL_INCLUDES)" LFLAGS="$(PARALLEL_LFLAGS)" LIBS="$(PARALLEL_LIBS)" CXX=$(PCC) ADMEXEC=$(PARALLEL_ADMEXEC) BAYESLIB_RULES=parallel BAYESLIB_NAME=parabayeslib.a 

define_admixmap:
#	touch config.h
	@echo "#define __ADMIXMAP__" > config.h
define_hapmixmap:
#	touch config.h
	@echo "#define __HAPMIXMAP__" > config.h

admixmap_message:	
	@echo **Compiling ADMIXMAP**

hapmixmap_message:	
	@echo **Compiling HAPMIXMAP**

##checks for existence of bayeslib. It it does not, compile it; if it does, write config.h
##this is necessary in case the version of config.h that exists is the wrong one
checkbayeslib:
	@echo **Checking for bayeslib**
	if test -r $(BAYESLIB_PATH)/$(BAYESLIB_NAME);then echo $(PARALLEL_DEFINES) >config.h;else $(MAKE) -C $(BAYESLIB_PATH) $(BAYESLIB_RULES) -e -f manual.make;fi;

bayeslib:
	$(MAKE) -C $(BAYESLIB_PATH) -e -fmanual.make

%.o: %.cc#rule for compilation
	$(CXX) $(CPPFLAGS) -I$(BAYESLIB_PATH) $(INCLUDES) -c $< -o $@

#recompile from scratch (serial)
new: 		clean all

#recompile everything from scratch (serial)
allnew:	libnew new

#delete all object files
clean:		          
	rm -f *.o config.h

#delete bayeslib object files
libclean: 
	$(MAKE) -C $(BAYESLIB_PATH) -fmanual.make clean

#delete sampler library
librealclean:
	$(MAKE) -C $(BAYESLIB_PATH) -fmanual.make realclean

#recompile sampler library
libnew: 
	@echo **recompiling bayeslib**
	$(MAKE) -C $(BAYESLIB_PATH) -e -fmanual.make allnew bayeslib

allclean: libclean clean

#delete all object files, libraries and exec
realclean: clean librealclean 
	rm -f $(DESTDIR)/$(EXEC)

#compile then run batch test script
check:		all
		perl $(DESTDIR)/batchtest.pl




