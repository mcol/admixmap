BASELIB_PATH=$(top_builddir)/src/base

AM_CPPFLAGS=-I$(BASELIB_PATH)

AM_LDFLAGS=-all-static -L$(top_builddir)/src/base

SUFFIXES=.cc .o
AM_CXXFLAGS=$(WFLAGS) -Wextra

bin_PROGRAMS = pedigree-tester

pedigree_tester_SOURCES = \
	pedigree-tester.cc


#pedigree_tester_LDADD=-lbase -lgsl -lgslcblas -lm
pedigree_tester_LDADD=-lbase -lm

pedigree_tester_DEPENDENCIES = $(top_builddir)/src/base/libbase.a



.PHONY: relink
relink:
	-$(RM) $(bin_PROGRAMS)
	$(MAKE)
