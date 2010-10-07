AM_CPPFLAGS = -I$(srcdir)/..
LDFLAGS=
SUFFIXES=.cc .o
OPTFLAGS=-fexceptions
AM_CXXFLAGS=$(WFLAGS) $(OPTFLAGS) -std=c++0x


if USE_OPENMP
    AM_CXXFLAGS+=-fopenmp
endif

#sources for the sublibrary
all_SOURCES = functions.cc DataReader.cc DataMatrix.cc LogWriter.cc OptionReader.cc StringConvertor.cc StringSplitter.cc GSLErrorHandler.cc ProbVector.cc pvector.cc IsNegative.h DelimitedBuf.cc Delimitedostream.cc RObjectWriter.cc exceptions.cc estr.cc

if PARALLEL
noinst_LTLIBRARIES = libutils-para.la
libutils_para_la_SOURCES=$(all_SOURCES)
libutils_para_la_CPPFLAGS="-DPARALLEL $(AM_CPPFLAGS)"

else
#this sublibrary will not be installed
noinst_LTLIBRARIES = libutils.la
libutils_la_SOURCES=$(all_SOURCES)
endif
#.cc.o:
#	$(CXX) $(AM_CXXFLAGS) $(CXXFLAGS) $(AM_CPPFLAGS) $(CPPFLAGS) -c $<



