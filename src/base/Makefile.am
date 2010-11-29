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
	AdmixtureProportions.h AdmixtureProportions.cc		\
	AlleleArray.h AlleleArray.cc				\
	AlleleFreqs.h AlleleFreqs.cc				\
	AlleleFreqParser.h AlleleFreqParser.cc			\
	AlleleFreqSampler.h AlleleFreqSampler.cc		\
	AncestryVector.h AncestryVector.cc			\
	Annealer.h Annealer.cc					\
	Chromosome.h Chromosome.cc				\
	CodeTimer.h CodeTimer.cc				\
	common.h common.cc					\
	CompositeLocus.h CompositeLocus.cc			\
	CopyNumberAssocTest.h CopyNumberAssocTest.cc		\
	DataValidError.h DataValidError.cc			\
	Filenames.h						\
	FreqArrays.h						\
	GeneticDistance.h GeneticDistance.cc			\
	GeneticDistanceUnit.h GeneticDistanceUnit.cc		\
	Genome.h Genome.cc					\
	Genotype.h Genotype.cc					\
	GenotypeIterator.h GenotypeIterator.cc			\
	GenotypeParser.h GenotypeParser.cc			\
	GenotypeProbOutputter.h GenotypeProbOutputter.cc	\
	GFileLexer.h GFileLexer.cc				\
	GPI.h GPI.cc						\
	Haplotype.h Haplotype.cc				\
	HapPair.h						\
	HiddenMarkovModel.h HiddenMarkovModel.cc		\
	HiddenMarkovModel.new.h HiddenMarkovModel.new.cc	\
	HiddenStateSpace.h HiddenStateSpace.cc			\
	HVIterator.h HVIterator.cc				\
	HWTest.h HWTest.cc					\
	Individual.h Individual.cc				\
	IndividualCollection.h IndividualCollection.cc		\
	InheritanceVector.h InheritanceVector.cc		\
	InputData.h InputData.cc				\
	misc.cc							\
	Model.h Model.cc					\
	Options.h Options.cc					\
	Organism.h Organism.cc					\
	OrganismArray.h						\
	PairSecIter.h						\
	PedBase.h PedBase.cc					\
	Pedigree.h Pedigree.cc					\
	PedigreeAddl.h PedigreeAddl.cc				\
	PedigreeGenStates.cc					\
	ProposalRingBuffer.h					\
	ResidualLDTest.h ResidualLDTest.cc			\
	RhoType.h						\
	ScoreTestBase.h ScoreTestBase.cc			\
	SimpleLocus.h SimpleLocus.cc				\
	SimpleLocusArray.h SimpleLocusArray.cc			\
	SimpleLocusParser.h SimpleLocusParser.cc		\
	TransProbCache.h TransProbCache.cc			\
	TwoDimArray.h						\
	$(NULL)
