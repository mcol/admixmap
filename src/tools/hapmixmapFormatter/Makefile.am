BAYESLIB_PATH = $(top_srcdir)/src/bcppcl
AM_CPPFLAGS=-I$(BAYESLIB_PATH)
AM_LDFLAGS = -L$(BAYESLIB_PATH)
noinst_PROGRAMS = FPHD FUD
FPHD_SOURCES = \
  FormatPhasedData.cpp \
  GenotypeEncoding.cpp \
  GenotypeEncoding.h \
  HapMapLegend.cpp \
  HapMapLegend.h \
  FPHDOptions.cpp \
  FPHDOptions.h \
  $(BAYESLIB_PATH)/utils/StringSplitter.cc

FUD_SOURCES = FormatUnphasedData.cpp

#TODO: only requires gsl_permutation.la sublibrary of gsl so wasteful to link against entire gsl and cblas
FPHD_LDADD=-lgsl -lgslcblas

EXTRA_DIST = repairSecondColumn.pl

#.cpp.o:
#$(CXX) $(AM_CXXFLAGS) $(AM_CPPFLAGS) -c $< -o $@
#.o:
#	$(CXX) $(AM_CXXFLAGS) $< -o $@

