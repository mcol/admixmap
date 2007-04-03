AM_CPPFLAGS = -I$(top_srcdir) -I$(top_srcdir)/samplers -I$(top_srcdir)/utils
LDFLAGS=
SUFFIXES=.cc .o
WFLAGS=-W -Wall
OPTFLAGS=
AM_CXXFLAGS=$(WFLAGS) $(OPTFLAGS)

noinst_LIBRARIES = libregression.a
libregression_a_SOURCES = Regression.cc LogisticRegression.cc LinearRegression.cc CoxRegression.cc
include_HEADERS = Regression.h LogisticRegression.h LinearRegression.h CoxRegression.h
includedir=$(prefix)/include/bcppcl

#.cc.o:
#	$(CXX) $(AM_CXXFLAGS) $(CXXFLAGS) $(AM_CPPFLAGS) $(CPPFLAGS) -c $<


