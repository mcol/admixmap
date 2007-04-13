BAYESLIB_PATH=$(srcdir)/../bcppcl
BASELIB_PATH=$(srcdir)/../base
AM_CPPFLAGS=-I$(BAYESLIB_PATH)  -I$(BASELIB_PATH)

AM_LDFLAGS= -L$(top_builddir)/src/bcppcl
SUFFIXES=.cc .o
WFLAGS=-W -Wall
AM_CXXFLAGS=$(WFLAGS)

ACLOCAL_AMFLAGS = -I m4 -I ../m4
#use -Wno-syntax if Automake is >= v1.10 to suppress "underquoted" warnings
#ACLOCAL_AMFLAGS = -Wno-syntax

all_SOURCES = HapMixFreqs.cc HapMixFreqs.h \
HapMixGenome.cc HapMixGenome.hh \
HapMixHMM.cc HapMixHMM.hh \
HapMixIndividual.cc HapMixIndividual.h \
HapMixIndividualCollection.cc HapMixIndividualCollection.h \
hapmixmap.cc \
HapMixModel.cc HapMixModel.h \
HapMixOptions.cc HapMixOptions.h \
MantelHaenszelTest.cc MantelHaenszelTest.h \
PopHapMix.cc PopHapMix.h 

parallel_SOURCES = 

if PARALLEL
bin_PROGRAMS = hapmixmap-para
hapmixmap_para_SOURCES = $(all_SOURCES) $(parallel_SOURCES)
hapmixmap_para_LDADD = $(top_builddir)/src/base/libbase-para.a -lbcppcl-para -lgsl -lgslcblas -lm $(PARALLEL_LIBS)
hapmixmap_para_CPPFLAGS = $(AM_CPPFLAGS) -DPARALLEL
hapmixmap_para_DEPENDENCIES=$(top_builddir)/src/base//libbase-para.a
else
bin_PROGRAMS = hapmixmap
hapmixmap_SOURCES = $(all_SOURCES)
hapmixmap_LDFLAGS = -static $(AM_LDFLAGS)
hapmixmap_LDADD = $(top_builddir)/src/base/libbase.a -lbcppcl -lgsl -lgslcblas -lm
hapmixmap_DEPENDENCIES=$(top_builddir)/src/base/libbase.a 
endif

#disable building of baselib with 'make distcheck'
DISTCHECK_CONFIGURE_FLAGS = --disable-baselib

#extra files to clean with 'make distclean'
#includes m4 files copied from elsewhere
DISTCLEANFILES = m4/*.m4

#extra files to distribute
#EXTRA_DIST = m4

#.cc.o:
#	$(CXX) $(AM_CXXFLAGS) $(CXXFLAGS) $(AM_CPPFLAGS) $(CPPFLAGS) -c $<

#bcppcl:
#	@(cd $(BAYESLIB_PATH) && $(MAKE) $(AM_MAKEFLAGS)) || exit 1;	
#base:
#	@(cd $(BASELIB_PATH) && $(MAKE) $(AM_MAKEFLAGS)) || exit 1;

