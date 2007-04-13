BAYESLIB_PATH=$(srcdir)/../bcppcl
BASELIB_PATH=$(srcdir)/../base
AM_CPPFLAGS=-I$(BAYESLIB_PATH)  -I$(BASELIB_PATH)

AM_LDFLAGS=-static -L$(top_builddir)/src/bcppcl
SUFFIXES=.cc .o
WFLAGS=-W -Wall
AM_CXXFLAGS=$(WFLAGS)

ACLOCAL_AMFLAGS = -I ../m4
#use -Wno-syntax if Automake is >= v1.10 to suppress "underquoted" warnings
#ACLOCAL_AMFLAGS = -Wno-syntax

bin_PROGRAMS = admixmap
admixmap_SOURCES = AdmixedIndividual.cc AdmixedIndividual.h \
AdmixIndividualCollection.cc AdmixIndividualCollection.h \
admixmap.cc \
AdmixMapModel.cc AdmixMapModel.h \
chib.cc chib.h \
DispersionFreqs.cc DispersionFreqs.h \
DispersionTest.cc DispersionTest.h \
IndAdmixOutputter.cc IndAdmixOutputter.h \
MisSpecAlleleFreqTest.cc MisSpecAlleleFreqTest.h \
PopAdmix.cc PopAdmix.h \
StratificationTest.cc StratificationTest.h

admixmap_LDADD=$(top_builddir)/src/base/libbase.a -lbcppcl -lgsl -lgslcblas -lm

admixmap_DEPENDENCIES = $(top_builddir)/src/base/libbase.a

#disable building of baselib with 'make distcheck'
DISTCHECK_CONFIGURE_FLAGS = --disable-baselib

#extra files to clean with 'make distclean'
#includes m4 files copied from elsewhere
DISTCLEANFILES = m4/*.m4

#.cc.o:
#	$(CXX) $(AM_CXXFLAGS) $(CXXFLAGS) $(AM_CPPFLAGS) $(CPPFLAGS) -c $<

#bcppcl:
#	@(cd $(BAYESLIB_PATH) && $(MAKE) $(AM_MAKEFLAGS)) || exit 1;	
#baselib:
#	@(cd $(BASELIB_PATH) && $(MAKE) $(AM_MAKEFLAGS)) || exit 1;

#$(BASELIB_PATH)/libbase.a : baselib

#lbbcppcl.a: bcppcl

allclean: clean
	@(cd $(BASELIB_PATH) && $(MAKE) $(AM_MAKEFLAGS) clean) || exit 1;

