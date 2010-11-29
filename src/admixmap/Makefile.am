BCLIB_PATH=$(srcdir)/../bclib
BASELIB_PATH=$(srcdir)/../base
AM_CPPFLAGS=-I$(BCLIB_PATH)  -I$(BASELIB_PATH)

CLEANFILES=svn_version.cc

AM_LDFLAGS=-all-static -L../bclib -L../base
SUFFIXES=.cc .o
AM_CXXFLAGS=$(WFLAGS) -Wextra

if USE_OPENMP
    AM_CXXFLAGS+=-fopenmp
endif

# This was needed in AM_CXXFLAGS when using reinterpret_cast<> for punning
# parent-IDs (rather than union):
#-Wstrict-aliasing=0


bin_PROGRAMS = admixmap
admixmap_SOURCES = \
	AdmixedIndividual.h AdmixedIndividual.cc			\
	AdmixFilenames.h						\
	AdmixFreqs.h AdmixFreqs.cc					\
	AdmixIndividualCollection.h AdmixIndividualCollection.cc	\
	AdmixMapModel.h AdmixMapModel.cc				\
	admixmap.cc							\
	AdmixmapGenotypeConverter.h AdmixmapGenotypeConverter.cc	\
	AdmixOptions.h AdmixOptions.cc					\
	AdmixPedigree.cc						\
	AdmixtureAssocTest.h AdmixtureAssocTest.cc			\
	AffectedsOnlyTest.h AffectedsOnlyTest.cc			\
	AllelicAssocSubTest.h AllelicAssocSubTest.cc			\
	AllelicAssocTest.h AllelicAssocTest.cc				\
	AncestryAssocTest.h AncestryAssocTest.cc			\
	chib.h chib.cc							\
	CorrelatedFreqs.h CorrelatedFreqs.cc				\
	DispersionFreqs.h DispersionFreqs.cc				\
	DispersionTest.h DispersionTest.cc				\
	IndAdmixOutputter.h IndAdmixOutputter.cc			\
	InputAdmixData.h InputAdmixData.cc				\
	MisSpecAlleleFreqTest.h MisSpecAlleleFreqTest.cc		\
	PopAdmix.h PopAdmix.cc						\
	ScoreTests.h ScoreTests.cc					\
	StratificationTest.h StratificationTest.cc			\
	$(NULL)

admixmap_LDADD=svn_version.o
admixmap_DEPENDENCIES=svn_version.o

admixmap_LDADD+=-lbase -lbcl -lgsl -lgslcblas -lm

if USE_OPENMP
    admixmap_LDADD+=-lgomp
endif


admixmap_DEPENDENCIES += $(top_builddir)/$(subdir)/../base/libbase.a
#admixmap_DEPENDENCIES += $(top_builddir)/src/base/libbase.a $(top_builddir)/src/bclib/libbcl.la

#.cc.o:
#	$(CXX) $(AM_CXXFLAGS) $(CXXFLAGS) $(AM_CPPFLAGS) $(CPPFLAGS) -c $<



svn_version.o: $(admixmap_SOURCES) Makefile.am
	svnversion $(srcdir) | sed 's/^/extern const char SVN_VERSION [] = "/;s/$$/";/' > $(srcdir)/svn_version.cc
	$(CXX) $(CXXFLAGS) -o "$@" -c $(srcdir)/svn_version.cc
