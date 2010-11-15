AM_CPPFLAGS = -I$(srcdir)/../bclib
AM_LDFLAGS= -L$(srcdir)/../bclib
SUFFIXES=.cc .o
#AM_CXXFLAGS=$(WFLAGS) -Wextra -DAGGRESSIVE_RANGE_CHECK=$(AGGRESSIVE_RANGE_CHECK)
AM_CXXFLAGS=$(WFLAGS) -Wextra

if USE_OPENMP
    AM_CXXFLAGS+=-fopenmp
endif


# This was needed in AM_CXXFLAGS when using reinterpret_cast<> for punning
# parent-IDs (rather than union):
#-Wstrict-aliasing=0


noinst_LIBRARIES = libbase.a

libbase_a_SOURCES = \
HVIterator.cc HVIterator.h \
AdmixtureProportions.cc AdmixtureProportions.h \
AlleleFreqs.cc AlleleFreqs.h \
AlleleFreqSampler.cc AlleleFreqSampler.h \
Annealer.cc Annealer.h \
Chromosome.cc Chromosome.h \
common.cc common.h \
CompositeLocus.cc CompositeLocus.h \
CopyNumberAssocTest.cc CopyNumberAssocTest.h \
Filenames.h \
FreqArrays.h \
GeneticDistanceUnit.cc GeneticDistanceUnit.h \
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
Genotype.h Genotype.cc \
GenotypeIterator.cc \
RhoType.h


if USE_GENOTYPE_PARSER
  libbase_a_SOURCES += \
		DataValidError.cc DataValidError.h \
		GFileLexer.cc GFileLexer.h \
		GeneticDistance.h GeneticDistance.cc \
		SimpleLocus.h SimpleLocus.cc \
		SimpleLocusArray.h SimpleLocusArray.cc \
		SimpleLocusParser.h SimpleLocusParser.cc \
		Organism.h Organism.cc \
		GenotypeParser.cc GenotypeParser.h \
		InheritanceVector.h InheritanceVector.cc \
		Pedigree.h Pedigree.cc PedigreeGenStates.cc \
		PedigreeAddl.h PedigreeAddl.cc \
		AncestryVector.cc AncestryVector.h \
		HiddenStateSpace.h HiddenStateSpace.cc \
		AlleleArray.h AlleleArray.cc \
		AlleleFreqParser.h AlleleFreqParser.cc \
		TransProbCache.h TransProbCache.cc \
		HiddenMarkovModel.new.h HiddenMarkovModel.new.cc \
		OrganismArray.h \
		PairSecIter.h \
		ProposalRingBuffer.h \
		TwoDimArray.h \
		PedBase.h PedBase.cc \
		CodeTimer.h CodeTimer.cc
else
  libbase_a_SOURCES += \
		GenotypeLoader.cc GenotypeLoader.h
endif


#.cc.o:
#	$(CXX) $(AM_CXXFLAGS) $(CXXFLAGS) $(AM_CPPFLAGS) $(CPPFLAGS) -c $<
