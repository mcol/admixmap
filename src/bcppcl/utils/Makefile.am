AM_CPPFLAGS = -I$(top_srcdir)
LDFLAGS=
SUFFIXES=.cc .o
WFLAGS=-W -Wall -fexceptions
OPTFLAGS=
AM_CXXFLAGS=$(WFLAGS) $(OPTFLAGS)

#this sublibrary will not be installed
noinst_LIBRARIES = libutils.a

#sources for the sublibrary
libutils_a_SOURCES = functions.cc DataReader.cc DataMatrix.cc LogWriter.cc StringConvertor.cc StringSplitter.cc GSLErrorHandler.cc  ProbVector.cc pvector.cc RObjectWriter.cc

#headers to be installed
include_HEADERS = ColumnIter.h DataMatrix.h DataReader.h dist.h Exceptions.h GSLExceptions.h GSLErrorHandler.h linalg.h misc.h LogWriter.h Array.hh ProbVector.h pvector.h  RObjectWriter.h SnapToZero.h StringSplitter.h StringConvertor.h

#place headers in bcppcl dir
includedir=$(prefix)/include/bcppcl

#.cc.o:
#	$(CXX) $(AM_CXXFLAGS) $(CXXFLAGS) $(AM_CPPFLAGS) $(CPPFLAGS) -c $<



