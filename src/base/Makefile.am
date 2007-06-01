AM_CPPFLAGS = -I$(srcdir)/../bcppcl
AM_LDFLAGS= -L$(srcdir)/../bcppcl
SUFFIXES=.cc .o
AM_CXXFLAGS=$(WFLAGS)

noinst_LIBRARIES = libbase.a

libbase_a_SOURCES = AlleleFreqs.cc AlleleFreqs.h \
AlleleFreqSampler.cc AlleleFreqSampler.h \
Annealer.cc Annealer.h \
Chromosome.cc Chromosome.h \
common.h \
Comms.cc Comms.h \
CompositeLocus.cc CompositeLocus.h \
CopyNumberAssocTest.cc CopyNumberAssocTest.h \
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
GenotypeLoader.cc GenotypeLoader.h 

#.cc.o:
#	$(CXX) $(AM_CXXFLAGS) $(CXXFLAGS) $(AM_CPPFLAGS) $(CPPFLAGS) -c $<


