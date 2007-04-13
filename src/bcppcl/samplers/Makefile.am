AM_CPPFLAGS = -I$(srcdir)/..
LDFLAGS=
SUFFIXES=.cc .o
WFLAGS=-W -Wall
OPTFLAGS=
AM_CXXFLAGS=$(WFLAGS) $(OPTFLAGS)

all_SOURCES = AdaptiveRejection.cc DirichletParamSampler.cc DispersionSampler.cc Gaussian.cc GaussianProposalMH.cc HamiltonianMonteCarlo.cc ModeFinder.cc MuSampler.cc rand.cc StepSizeTuner.cc 

if PARALLEL
noinst_LTLIBRARIES = libsampler-para.la
libsampler_para_la_SOURCES	= $(all_SOURCES)
libsampler_para_la_CPPFLAGS	= $(AM_CPPFLAGS) -DPARALLEL
else
noinst_LTLIBRARIES = libsampler.la
libsampler_la_SOURCES = $(all_SOURCES)
endif
#.cc.o:
#	$(CXX) $(AM_CXXFLAGS) $(CXXFLAGS) $(AM_CPPFLAGS) $(CPPFLAGS) -c $<


