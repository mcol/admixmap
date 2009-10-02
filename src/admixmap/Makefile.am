BCLIB_PATH=$(srcdir)/../bclib
BASELIB_PATH=$(srcdir)/../base
AM_CPPFLAGS=-I$(BCLIB_PATH)  -I$(BASELIB_PATH)

CLEANFILES=svn_version.cc

AM_LDFLAGS=-all-static -L$(top_builddir)/src/bclib -L$(top_builddir)/src/base
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
AdmixPedigree.cc \
AdmixFilenames.h \
AdmixedIndividual.cc AdmixedIndividual.h \
AdmixIndividualCollection.cc AdmixIndividualCollection.h \
admixmap.cc \
AdmixMapModel.cc AdmixMapModel.h \
AdmixOptions.cc AdmixOptions.h \
AdmixFreqs.cc AdmixFreqs.h \
chib.cc chib.h \
CorrelatedFreqs.cc CorrelatedFreqs.h \
DispersionFreqs.cc DispersionFreqs.h \
DispersionTest.cc DispersionTest.h \
IndAdmixOutputter.cc IndAdmixOutputter.h \
InputAdmixData.cc InputAdmixData.h \
MisSpecAlleleFreqTest.cc MisSpecAlleleFreqTest.h \
PopAdmix.cc PopAdmix.h \
StratificationTest.cc StratificationTest.h \
AffectedsOnlyTest.cc AffectedsOnlyTest.h \
AncestryAssocTest.cc AncestryAssocTest.h \
AdmixtureAssocTest.cc AdmixtureAssocTest.h \
AllelicAssocTest.cc AllelicAssocTest.h \
AllelicAssocSubTest.cc AllelicAssocSubTest.h \
ScoreTests.cc ScoreTests.h

if USE_GENOTYPE_PARSER
  admixmap_SOURCES += \
    AdmixmapGenotypeConverter.h AdmixmapGenotypeConverter.cc
endif


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
	svnversion $(top_builddir)/src | sed 's/^/extern const char SVN_VERSION [] = "/;s/$$/";/' > $(srcdir)/svn_version.cc
	$(CXX) $(CXXFLAGS) -o "$@" -c $(srcdir)/svn_version.cc
