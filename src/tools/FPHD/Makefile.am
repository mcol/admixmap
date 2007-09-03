AM_CPPFLAGS=-I$(top_srcdir)/src/bclib
AM_LDFLAGS = -all-static -L$(top_builddir)/src/bclib

bin_PROGRAMS = FPHD
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

#TODO: only requires gsl_permutation.la sublibrary of gsl so wasteful to link against entire gsl
FPHD_LDADD= -lbcl -lgsl

EXTRA_DIST = getdata.pl README.txt

