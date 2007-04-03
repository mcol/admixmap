SUBDIRS=utils samplers regression
#SUBLIBS=utils samplers regression
SUBLIBS=utils/*.o samplers/*.o regression/*.o

#name of the library
lib_LIBRARIES=libbcppcl.a

#add all sublibraries
libbcppcl_a_LIBADD=$(SUBLIBS)

#no sources in this directory. Must be defined as blank or automake will look for bcppcl.c
libbcppcl_a_SOURCES=
.SUFFIXES=

#headers to be installed
include_HEADERS = bcppcl.h common.h

#install to this directory
includedir=$(prefix)/include/bcppcl

#local m4 macros and flags for aclocal
#use -Wno-syntax to suppress 'underquoted definition' warnings
#but that's only valid for aclocal >= v1.10
ACLOCAL_AMFLAGS = -I m4
#ACLOCAL_AMFLAGS = -I m4-Wno-syntax

#extra files to include
EXTRA_DIST = INSTALL

#build sampler library
bsl:
	@(cd samplers && $(MAKE) $(AM_MAKEFLAGS)) || exit 1;

#build regression library
brl:
	@(cd regression && $(MAKE) $(AM_MAKEFLAGS)) || exit 1;

#build utility library
bul:
	@(cd utils && $(MAKE) $(AM_MAKEFLAGS)) || exit 1;