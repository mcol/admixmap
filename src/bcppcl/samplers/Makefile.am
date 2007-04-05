AM_CPPFLAGS = -I$(top_srcdir)
LDFLAGS=
SUFFIXES=.cc .o
WFLAGS=-W -Wall
OPTFLAGS=
AM_CXXFLAGS=$(WFLAGS) $(OPTFLAGS)

noinst_LTLIBRARIES = libsampler.la
libsampler_la_SOURCES = AdaptiveRejection.cc DirichletParamSampler.cc DispersionSampler.cc Gaussian.cc GaussianProposalMH.cc HamiltonianMonteCarlo.cc ModeFinder.cc MuSampler.cc rand.cc StepSizeTuner.cc 

include_HEADERS = AdaptiveRejection.h DirichletParamSampler.h DispersionSampler.h Gaussian.h GaussianProposalMH.h HamiltonianMonteCarlo.h ModeFinder.h MuSampler.h rand.h gsl-sprng.h StepSizeTuner.h
includedir=$(prefix)/include/bcppcl

#.cc.o:
#	$(CXX) $(AM_CXXFLAGS) $(CXXFLAGS) $(AM_CPPFLAGS) $(CPPFLAGS) -c $<


