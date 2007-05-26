BAYESLIB_PATH=$(srcdir)/../bcppcl
BASELIB_PATH=$(srcdir)/../base
AM_CPPFLAGS=-I$(BAYESLIB_PATH)  -I$(BASELIB_PATH)

AM_LDFLAGS= -L$(top_builddir)/src/bcppcl
SUFFIXES=.cc .o
AM_CXXFLAGS=$(WFLAGS)

all_SOURCES = hapmixmap.cc \
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

parallel_SOURCES = 

if PARALLEL
bin_PROGRAMS = hapmixmap-para
hapmixmap_para_SOURCES = $(all_SOURCES) $(parallel_SOURCES)
hapmixmap_para_LDADD = $(top_builddir)/$(subdir)/../base/libbase-para.a -lbcppcl-para -lgsl -lgslcblas -lm $(PARALLEL_LIBS)
hapmixmap_para_CPPFLAGS = $(AM_CPPFLAGS) -DPARALLEL
hapmixmap_para_DEPENDENCIES=$(top_builddir)/$(subdir)/../base/libbase-para.a
else
bin_PROGRAMS = hapmixmap
hapmixmap_SOURCES = $(all_SOURCES)
hapmixmap_LDFLAGS = -all-static $(AM_LDFLAGS)
hapmixmap_LDADD = $(top_builddir)/$(subdir)/../base/libbase.a -lbcppcl -lgsl -lgslcblas -lm
hapmixmap_DEPENDENCIES=$(top_builddir)/$(subdir)/../base/libbase.a 
endif

#.cc.o:
#	$(CXX) $(AM_CXXFLAGS) $(CXXFLAGS) $(AM_CPPFLAGS) $(CPPFLAGS) -c $<

