noinst_LTLIBRARIES         	= libutils.la
libutils_la_SOURCES		= functions.cc \
				  DataReader.cc \
				  DataMatrix.cc \
				  StringSplitter.cc \
				  StringConvertor.cc \
				  LogWriter.cc \
				  RObjectWriter.cc \
				  GSLErrorHandler.cc \
				  ProbVector.cc \
			          pvector.cc \
				  SnapToZero.cc \
				  DivideBy.cc \
				  ZeroOrMore.cc
>>>>>>> .r2243

noinst_HEADERS			= DataMatrix.h \
				  DataReader.h \
				  dist.h \
				  linalg.h \
				  LogWriter.h \
				  misc.h \
				  StringConvertor.h \
				  StringSplitter.h \
				  ColumnIter.h \
				  RObjectWriter.h
				  GSLErrorHandler.h \
				  Exceptions.h \
				  GSLExceptions.h \
                                  ProbVector.h \
				  pvector.h \
				  SnapToZero.h \
				  DivideBy.h \
				  ZeroOrMore.h

INCLUDES			= -I$(top_srcdir)/src/bayeslib

if PARALLEL
noinst_LTLIBRARIES		+= libutils-para.la
libutils_para_la_SOURCES	= $(libutils_la_SOURCES)
libutils_para_la_CPPFLAGS	= -DPARALLEL
endif

# libutils_la_LDFLAGS		= -version-info 0:0:0
