AUTOMAKE_OPTIONS	= gnu
LIBS			= -lgsl -lgslcblas $(PARALLEL_LIBS)
# LDADD			= libgenepi.la

INCLUDES		= -I$(top_srcdir)/src/bayeslib


bin_PROGRAMS		= hapmixmap admixmap
noinst_LTLIBRARIES	= libgenepi.la

hapmixmap_LDFLAGS	= $(STATIC_LDFLAGS)
hapmixmap_CPPFLAGS	= -D__HAPMIXMAP__
admixmap_CPPFLAGS	= -D__ADMIXMAP__

hapmixmap_SOURCES	= hapmixmap.cc \
			HapMixFreqs.cc \
		   	Model.cc \
			HapMixIndividual.cc \
			HapMixIndividualCollection.cc \
			HapMixModel.cc \
			HapMixOptions.cc \
			MantelHaenszelTest.cc \
			PopHapMix.cc

hapmixmap_LDADD		= libgenepi.la

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

admixmap_LDADD		= libgenepi.la


noinst_HEADERS	= admixmap.h \
		  HapMixFreqs.h \
			  HapMixIndividual.h \
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
			  HapPair.h \
			  Haplotype.h \
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
			  ScoreTests.h \
			  GPI.h \
			  GenotypeIterator.h \
			  FreqArrays.h


libgenepi_la_SOURCES = misc.cc \
		       AlleleFreqs.cc \
		       Options.cc \
		       AdmixOptions.cc \
		       InputData.cc \
		       Genome.cc \
		       Chromosome.cc \
		       Haplotype.cc \
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

libgenepi_la_LIBADD = $(top_builddir)/src/bayeslib/libbayes.la

if PARALLEL

# Parallel libgenepi
noinst_LTLIBRARIES		+= libgenepi-para.la
libgenepi_para_la_SOURCES	= $(libgenepi_la_SOURCES)
libgenepi_para_la_CPPFLAGS	= -DPARALLEL
libgenepi_para_la_LIBADD	= $(top_builddir)/src/bayeslib/libbayes-para.la

# Parallel hapmixmap
bin_PROGRAMS			+= hapmixmap-para
hapmixmap_para_CPPFLAGS		= $(hapmixmap_CPPFLAGS) -DPARALLEL
hapmixmap_para_SOURCES		= $(hapmixmap_SOURCES)
# hapmixmap_para_LIBS		= $(PARALLEL_LIBS)
LIBS				+= $(PARALLEL_LIBS)
hapmixmap_para_LDADD		= $(top_builddir)/src/bayeslib/libbayes-para.la \
				  libgenepi-para.la
endif
