# This works fine when compiling "in-place" (in the source directory) but not
# when building releases in a separate directory:
BASELIB_OBJ_PATH=$(top_builddir)/src/base
BASELIB_SRC_PATH=$(srcdir)/../../../src/base

AM_CPPFLAGS=-I$(BASELIB_SRC_PATH)

AM_LDFLAGS=-all-static -L$(BASELIB_OBJ_PATH)

SUFFIXES=.cc .o
AM_CXXFLAGS=$(WFLAGS) -Wextra

if USE_OPENMP
    AM_CXXFLAGS+=-fopenmp
endif

bin_PROGRAMS = pedigree-tester

pedigree_tester_SOURCES = \
	pedigree-tester.cc


#pedigree_tester_LDADD=-lbase -lgsl -lgslcblas -lm
pedigree_tester_LDADD=-lbase -lm

if USE_OPENMP
    pedigree_tester_LDADD+=-lgomp
endif



pedigree_tester_DEPENDENCIES = $(BASELIB_OBJ_PATH)/libbase.a



.PHONY: relink
relink:
	-$(RM) $(bin_PROGRAMS)
	$(MAKE)
