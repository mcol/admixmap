BAYESLIB_PATH=$(srcdir)/../bcppcl
BASELIB_PATH=$(srcdir)/../base
AM_CPPFLAGS=-I$(BAYESLIB_PATH)  -I$(BASELIB_PATH)

AM_LDFLAGS=-static -L$(top_builddir)/src/bcppcl
SUFFIXES=.cc .o
WFLAGS=-W -Wall
AM_CXXFLAGS=$(WFLAGS)

bin_PROGRAMS = admixmap
admixmap_SOURCES = AdmixedIndividual.cc AdmixedIndividual.h \
AdmixIndividualCollection.cc AdmixIndividualCollection.h \
admixmap.cc \
AdmixMapModel.cc AdmixMapModel.h \
AdmixOptions.cc AdmixOptions.h \
chib.cc chib.h \
DispersionFreqs.cc DispersionFreqs.h \
DispersionTest.cc DispersionTest.h \
IndAdmixOutputter.cc IndAdmixOutputter.h \
MisSpecAlleleFreqTest.cc MisSpecAlleleFreqTest.h \
PopAdmix.cc PopAdmix.h \
StratificationTest.cc StratificationTest.h

admixmap_LDADD=$(builddir)/../base/libbase.a -lbcppcl -lgsl -lgslcblas -lm

admixmap_DEPENDENCIES = $(builddir)/../base/libbase.a

#.cc.o:
#	$(CXX) $(AM_CXXFLAGS) $(CXXFLAGS) $(AM_CPPFLAGS) $(CPPFLAGS) -c $<

