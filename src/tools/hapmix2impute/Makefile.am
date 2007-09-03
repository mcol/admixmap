FPHD_DIR = $(top_srcdir)/src/tools/FPHD
AM_CPPFLAGS=-I$(top_srcdir)/src/bclib -I$(FPHD_DIR)
AM_LDFLAGS = -all-static -L$(top_builddir)/src/bclib 
bin_PROGRAMS = hapmix2impute
hapmix2impute_SOURCES = hapmix2impute.cpp 


hapmix2impute_LDADD= -lbcl ../FPHD/HapMapLegend.o
EXTRA_DIST = hapmix2impute.R
