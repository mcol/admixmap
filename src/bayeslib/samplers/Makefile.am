noinst_LTLIBRARIES         	= libsampler.la

libsampler_la_SOURCES		= rand.cc \
				  AdaptiveRejection.cc \
				  StepSizeTuner.cc \
				  Gaussian.cc \
				  GaussianProposalMH.cc \
				  DirichletParamSampler.cc \
				  HamiltonianMonteCarlo.cc \
				  MuSampler.cc \
				  DispersionSampler.cc \
				  ModeFinder.cc

noinst_HEADERS		= AdaptiveRejection.h \
				  DirichletParamSampler.h \
				  DispersionSampler.h \
				  Gaussian.h \
				  GaussianProposalMH.h \
				  gsl-sprng.h \
				  HamiltonianMonteCarlo.h \
				  ModeFinder.h \
				  MuSampler.h \
				  rand.h \
				  StepSizeTuner.h
INCLUDES			= -I$(top_srcdir)/src/bayeslib

if PARALLEL
noinst_LTLIBRARIES		+= libsampler-para.la
libsampler_para_la_SOURCES	= $(libsampler_la_SOURCES)
libsampler_para_la_CPPFLAGS	= -DPARALLEL
endif

# libsampler_a_LDFLAGS		= -version-info 0:0:0
