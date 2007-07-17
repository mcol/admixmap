AM_CPPFLAGS = -I$(srcdir)/..
LDFLAGS=
SUFFIXES=.cc .o
OPTFLAGS=
AM_CXXFLAGS=$(WFLAGS) $(OPTFLAGS)

all_SOURCES = Regression.cc LogisticRegression.cc LinearRegression.cc CoxRegression.cc

if PARALLEL
noinst_LTLIBRARIES = libregression-para.la
libregression_para_la_SOURCES=$(all_SOURCES)

else
noinst_LTLIBRARIES = libregression.la
libregression_la_SOURCES=$(all_SOURCES)
endif

#.cc.o:
#	$(CXX) $(AM_CXXFLAGS) $(CXXFLAGS) $(AM_CPPFLAGS) $(CPPFLAGS) -c $<


