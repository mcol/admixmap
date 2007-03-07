BAYESLIB_PATH=../../bayeslib
AM_CPPFLAGS=$(BAYESLIB_PATH)/utils
bin_PROGRAMS = DataFormatter
SUFFIXES=.cc .cpp .o
DataFormatter_SOURCES = FormatPhasedData.cpp GenotypeEncoder.cpp $(BAYESLIB_PATH)/utils/StringSplitter.cc

#TODO: only requires gsl_permutation.la sublibrary of gsl so wasteful to link against entire gsl and cblas
LDADD=-lgsl -lgslcblas
