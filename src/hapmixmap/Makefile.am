BAYESLIB_PATH=$(srcdir)/../bcppcl
BASELIB_PATH=$(srcdir)/../base
AM_CPPFLAGS=-I$(BAYESLIB_PATH)  -I$(BASELIB_PATH)

AM_LDFLAGS= -L$(top_builddir)/src/bcppcl
SUFFIXES=.cc .o
AM_CXXFLAGS=$(WFLAGS) -fexceptions

bin_PROGRAMS = hapmixmap
hapmixmap_LDFLAGS = -all-static $(AM_LDFLAGS)
hapmixmap_LDADD = $(top_builddir)/$(subdir)/../base/libbase.a -lbcppcl -lgsl -lgslcblas -lm
hapmixmap_DEPENDENCIES=$(top_builddir)/$(subdir)/../base/libbase.a 

hapmixmap_SOURCES = hapmixmap.cc \
HapMixFilenames.h \
HapMixFreqs.cc HapMixFreqs.h \
HapMixGenome.cc HapMixGenome.hh \
HapMixHMM.cc HapMixHMM.hh \
HapMixIndividual.cc HapMixIndividual.h \
HapMixIndividualCollection.cc HapMixIndividualCollection.h \
HapMixModel.cc HapMixModel.h \
HapMixOptions.cc HapMixOptions.h \
HapMixAllelicAssocTest.cc HapMixAllelicAssocTest.h  \
HapMixGenotypeLoader.c HapMixGenotypeLoader.h \
InputHapMixData.c InputHapMixData.h \
MantelHaenszelTest.cc MantelHaenszelTest.h \
PopHapMix.cc PopHapMix.h

#.cc.o:
#	$(CXX) $(AM_CXXFLAGS) $(CXXFLAGS) $(AM_CPPFLAGS) $(CPPFLAGS) -c $<

