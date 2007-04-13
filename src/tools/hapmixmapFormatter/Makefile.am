BAYESLIB_PATH = $(top_srcdir)/src/bcppcl
AM_CPPFLAGS=-I$(BAYESLIB_PATH)
AM_LDFLAGS = -L$(BAYESLIB_PATH)
noinst_PROGRAMS = DataFormatter
DataFormatter_SOURCES = \
  FormatPhasedData.cpp \
  GenotypeEncoder.cpp \
  Formatting.h \
  ranker.h \
  $(BAYESLIB_PATH)/utils/StringSplitter.cc

#TODO: only requires gsl_permutation.la sublibrary of gsl so wasteful to link against entire gsl and cblas
LDADD=-lgsl -lgslcblas

#.cpp.o:
#$(CXX) $(AM_CXXFLAGS) $(AM_CPPFLAGS) -c $< -o $@
#.o:
#	$(CXX) $(AM_CXXFLAGS) $< -o $@

