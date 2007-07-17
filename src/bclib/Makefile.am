SUBDIRS=utils samplers regression bclib doxygen . examples
SUBLIBS=utils/libutils.la samplers/libsampler.la regression/libregression.la
PARALLEL_SUBLIBS=utils/libutils-para.la samplers/libsampler-para.la regression/libregression-para.la

#SUBLIBS=utils/*.o samplers/*.o regression/*.o

if PARALLEL
#name of the parallel library
lib_LTLIBRARIES=libbcl-para.la
#include all parallel sublibs
libbcl_para_la_LIBADD=$(PARALLEL_SUBLIBS)
libbcl_para_la_CPPFLAGS = $(AM_CPPFLAGS) -DPARALLEL
#no sources in this directory
libbcl_para_la_SOURCES=
else

#name of the serial library
lib_LTLIBRARIES=libbcl.la

#add all sublibraries
libbcl_la_LIBADD=$(SUBLIBS)

#no sources in this directory. Must be defined as blank or automake will look for bcppcl.c
libbcl_la_SOURCES=
endif

.SUFFIXES=

#local m4 macros and flags for aclocal
ACLOCAL_AMFLAGS = -I m4

#use -Wno-syntax if Automake is >= v1.10 to suppress "underquoted" warnings
#ACLOCAL_AMFLAGS = -I m4 -Wno-syntax


#extra files to distribute
EXTRA_DIST = INSTALL README.rtf

#skip distcleancheck - it always fails
distuninstallcheck:
	@:

#build sampler library
bsl:
	@(cd samplers && $(MAKE) $(AM_MAKEFLAGS)) || exit 1;

#build regression library
brl:
	@(cd regression && $(MAKE) $(AM_MAKEFLAGS)) || exit 1;

#build utility library
bul:
	@(cd utils && $(MAKE) $(AM_MAKEFLAGS)) || exit 1;

#create docs
docs:
	@(cd doxygen && $(MAKE) $(AM_MAKEFLAGS)) || exit 1;
