AM_CPPFLAGS = -I$(srcdir)/..
LDFLAGS=
SUFFIXES=.cc .o
WFLAGS=-W -Wall -fexceptions
OPTFLAGS=
AM_CXXFLAGS=$(WFLAGS) $(OPTFLAGS)


#sources for the sublibrary
all_SOURCES = functions.cc DataReader.cc DataMatrix.cc LogWriter.cc StringConvertor.cc StringSplitter.cc GSLErrorHandler.cc  ProbVector.cc pvector.cc RObjectWriter.cc

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



