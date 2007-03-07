SUBDIRS				= mock

TESTS				= run-tests
check_PROGRAMS			= $(TESTS)

run_tests_SOURCES		= run-tests.cc
run_tests_LDADD			= libtestconverter.la \
				  libtestindividual.la \
				  ../bayeslib/libbayes.la
run_tests_LDFLAGS		= $(CPPUNIT_LIBS) -lgsl -lgslcblas -lmockpp
run_tests_CXXFLAGS		= $(CPPUNIT_CFLAGS)

check_LTLIBRARIES		= libtestconverter.la \
				  libtestindividual.la

libtestconverter_la_SOURCES	= UserDataLocusTest.cc \
				  LegendTest.cpp \
				  ../tools/UserDataLocus.cc \
				  ../tools/Legend.cc \
				  ../tools/UserDataIndividual.cc \
				  ../tools/UserData.cc \
				  LegendData.h

libtestindividual_la_SOURCES	= HapMixIndividualTest.cc \
				  CompositeLocusTest.cc

libtestindividual_la_LIBADD	= mock/libmock.la ../admixmap/libhapmixmap.la

noinst_HEADERS			= HapMixIndividualTest.h \
				  UserDataLocusTest.h \
				  LegendTest.h \
				  CompositeLocusTest.h

INCLUDES			= -I $(top_srcdir)/src/bayeslib
