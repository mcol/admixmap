BAYESLIB_PATH=$(srcdir)/../bcppcl
BASELIB_PATH=$(srcdir)/../base
AM_CPPFLAGS=-I$(BAYESLIB_PATH)  -I$(BASELIB_PATH)

AM_LDFLAGS=-all-static -L$(top_builddir)/src/bcppcl
SUFFIXES=.cc .o
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
InputAdmixData.cc InputAdmixData.h \
MisSpecAlleleFreqTest.cc MisSpecAlleleFreqTest.h \
PopAdmix.cc PopAdmix.h \
StratificationTest.cc StratificationTest.h \
AffectedsOnlyTest.cc AffectedsOnlyTest.h \
AncestryAssocTest.cc AncestryAssocTest.h \
AdmixtureAssocTest.cc AdmixtureAssocTest.h \
AllelicAssocTest.cc AllelicAssocTest.h \
AllelicAssocSubTest.cc AllelicAssocSubTest.h \
ScoreTests.cc ScoreTests.h 

admixmap_LDADD=$(top_builddir)/$(subdir)/../base/libbase.a -lbcppcl -lgsl -lgslcblas -lm

admixmap_DEPENDENCIES = $(top_builddir)/$(subdir)/../base/libbase.a

#.cc.o:
#	$(CXX) $(AM_CXXFLAGS) $(CXXFLAGS) $(AM_CPPFLAGS) $(CPPFLAGS) -c $<

