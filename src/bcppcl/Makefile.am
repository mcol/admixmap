SUBDIRS=utils samplers regression doxygen
SUBLIBS=utils/libutils.la samplers/libsampler.la regression/libregression.la
#SUBLIBS=utils/*.o samplers/*.o regression/*.o

#name of the library
lib_LTLIBRARIES=libbcppcl.la

#add all sublibraries
libbcppcl_la_LIBADD=$(SUBLIBS)

#no sources in this directory. Must be defined as blank or automake will look for bcppcl.c
libbcppcl_la_SOURCES=
.SUFFIXES=

#headers to be installed
include_HEADERS = bcppcl.h common.h

#install to this directory
includedir=$(prefix)/include/bcppcl

#local m4 macros and flags for aclocal
ACLOCAL_AMFLAGS = -I m4 -Wno-syntax

#extra files to distribute
EXTRA_DIST = INSTALL README.rtf

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