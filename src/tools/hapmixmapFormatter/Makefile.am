BAYESLIB_PATH = $(top_srcdir)/src/bayeslib
AM_CPPFLAGS=-I$(BAYESLIB_PATH)/utils
bin_PROGRAMS = DataFormatter
DataFormatter_SOURCES = FormatPhasedData.cpp		\
			GenotypeEncoder.cpp 		\
			$(BAYESLIB_PATH)/utils/StringSplitter.cc

noinst_HEADERS = Formatting.h		\
		 Formatting2.h		\
		 ranker.h

#TODO: only requires gsl_permutation.la sublibrary of gsl so wasteful to link against entire gsl and cblas
LDADD=-lgsl -lgslcblas
