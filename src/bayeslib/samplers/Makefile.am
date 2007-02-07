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
INCLUDES			= -I ..
# libsampler_a_LDFLAGS		= -version-info 0:0:0
