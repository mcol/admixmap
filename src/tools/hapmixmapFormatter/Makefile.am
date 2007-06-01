AM_CPPFLAGS=-I$(top_srcdir)/src/bcppcl
AM_LDFLAGS = -all-static -L$(top_builddir)/src/bcppcl
noinst_PROGRAMS = FPHD FUD
FPHD_SOURCES = \
  FormatPhasedData.cpp \
  GenotypeEncoding.cpp \
  GenotypeEncoding.h \
  HapMapLegend.cpp \
  HapMapLegend.h \
  FPHDOptions.cpp \
  FPHDOptions.h \
  HapMapData.h \
  HapMapData.cpp \
  UserGenotypes.h \
  UserGenotypes.cpp \
  InitialValues.h \
  InitialValues.cpp


FUD_SOURCES = FormatUnphasedData.cpp

#TODO: only requires gsl_permutation.la sublibrary of gsl so wasteful to link against entire gsl
FPHD_LDADD= -lbcppcl -lgsl

EXTRA_DIST = repairSecondColumn.pl getdata.pl README.txt

#.cpp.o:
#$(CXX) $(AM_CXXFLAGS) $(AM_CPPFLAGS) -c $< -o $@
#.o:
#	$(CXX) $(AM_CXXFLAGS) $< -o $@

