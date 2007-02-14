AUTOMAKE_OPTIONS	= gnu
LIBS			= -lgsl -lgslcblas $(PARALLEL_LIBS)
LDADD			= libgenepi.la

# The following line should be removed for the package to conform with
# the packaging standard, and to allow builds outside the source
# directory.
INCLUDES		= -I ../bayeslib


bin_PROGRAMS		= hapmixmap admixmap
noinst_LTLIBRARIES	= libgenepi.la

hapmixmap_LDFLAGS	= $(STATIC_LDFLAGS)
hapmixmap_CPPFLAGS	= -D__HAPMIXMAP__
admixmap_CPPFLAGS	= -D__ADMIXMAP__

hapmixmap_SOURCES	= hapmixmap.cc \
			HapMixFreqs.cc \
		   	Model.cc \
			HapMixIndividualCollection.cc \
			HapMixModel.cc \
			HapMixOptions.cc \
			MantelHaenszelTest.cc \
			PopHapMix.cc


hapmixmap_LDADD		= ../bayeslib/libbayes.la \
			  libgenepi.la

admixmap_SOURCES = admixmap.cc \
		   Model.cc \
		   AdmixedIndividual.cc \
		   AdmixIndividualCollection.cc \
		   AdmixMapModel.cc \
		   chib.cc \
		   DispersionFreqs.cc \
		   DispersionTest.cc \
		   IndAdmixOutputter.cc \
		   MisSpecAlleleFreqTest.cc \
		   PopAdmix.cc \
		   StratificationTest.cc

admixmap_LDADD		= ../bayeslib/libbayes.la \
			  libgenepi.la


noinst_HEADERS	= admixmap.h \
		  HapMixFreqs.h \
			  HapMixIndividualCollection.h \
			  HapMixModel.h \
			  HapMixOptions.h \
			  AdmixedIndividual.h \
			  AdmixIndividualCollection.h \
			  AdmixMapModel.h \
			  AdmixOptions.h \
			  AdmixtureAssocTest.h \
			  chib.h \
			  DispersionFreqs.h \
			  DispersionTest.h \
			  IndAdmixOutputter.h \
			  MisSpecAlleleFreqTest.h \
			  PopAdmix.h \
			  StratificationTest.h \
			  AffectedsOnlyTest.h \
			  AlleleFreqSampler.h \
			  AlleleFreqs.h \
			  AncestryAssocTest.h \
			  Annealer.h \
			  Chromosome.h \
			  Comms.h \
			  CompositeLocus.h \
			  Genome.h \
			  GenotypeProbOutputter.h \
			  HMM.h \
			  HWTest.h \
			  IndividualCollection.h \
			  Individual.h \
			  InputData.h \
			  MantelHaenszelTest.h \
			  Model.h \
			  Options.h \
			  PopHapMix.h \
			  ResidualLDTest.h \
			  ScoreTestBase.h \
			  ScoreTests.h


libgenepi_la_SOURCES = misc.cc \
		       AlleleFreqs.cc \
		       Options.cc \
		       AdmixOptions.cc \
		       InputData.cc \
		       Genome.cc \
		       Chromosome.cc \
		       CompositeLocus.cc \
		       IndividualCollection.cc \
		       Individual.cc \
		       HMM.cc \
		       AdmixtureAssocTest.cc \
		       AffectedsOnlyTest.cc \
		       AncestryAssocTest.cc \
		       ScoreTestBase.cc \
		       ScoreTests.cc \
		       ResidualLDTest.cc \
		       HWTest.cc \
		       AlleleFreqSampler.cc \
		       Comms.cc \
		       GenotypeProbOutputter.cc \
		       Annealer.cc

libgenepi_la_LIBADD = ../bayeslib/libbayes.la


