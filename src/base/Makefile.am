AM_CPPFLAGS = -I$(srcdir)/../bcppcl
AM_LDFLAGS= -L$(srcdir)/../bcppcl
SUFFIXES=.cc .o
AM_CXXFLAGS=$(WFLAGS)

all_SOURCES = AdmixtureAssocTest.cc AdmixtureAssocTest.h \
AffectedsOnlyTest.cc AffectedsOnlyTest.h \
AlleleFreqs.cc AlleleFreqs.h \
AlleleFreqSampler.cc AlleleFreqSampler.h \
AllelicAssocTest.cc AllelicAssocTest.h \
AllelicAssocSubTest.cc AllelicAssocSubTest.h \
Annealer.cc Annealer.h \
Chromosome.cc Chromosome.h \
common.h \
Comms.cc Comms.h \
CompositeLocus.cc CompositeLocus.h \
CopyNumberAssocTest.cc CopyNumberAssocTest.h \
EventLogger.cc EventLogger.hh \
FreqArrays.h \
GeneticDistanceUnit.h \
Genome.cc Genome.h \
GenotypeIterator.h \
GenotypeProbOutputter.cc GenotypeProbOutputter.h \
GPI.cc GPI.h \
Haplotype.cc Haplotype.h \
HapPair.h \
HiddenMarkovModel.cc HiddenMarkovModel.h \
HWTest.cc HWTest.h \
Individual.cc Individual.h \
IndividualCollection.cc IndividualCollection.h \
InputData.cc InputData.h  \
misc.cc \
Model.cc Model.h \
Options.cc Options.h \
ResidualLDTest.cc ResidualLDTest.h \
ScoreTestBase.cc ScoreTestBase.h \
ScoreTests.cc ScoreTests.h 

if PARALLEL
noinst_LIBRARIES = libbase-para.a
libbase_para_a_SOURCES = $(all_SOURCES)
libbase_para_a_CPPFLAGS = $(AM_CPPFLAGS) -DPARALLEL
else
noinst_LIBRARIES = libbase.a
libbase_a_SOURCES = $(all_SOURCES)
endif


#.cc.o:
#	$(CXX) $(AM_CXXFLAGS) $(CXXFLAGS) $(AM_CPPFLAGS) $(CPPFLAGS) -c $<


