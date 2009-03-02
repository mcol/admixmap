
SUBDIRS = bclib base

if ! PARALLEL
  SUBDIRS += admixmap
endif

if ! USE_GENOTYPE_PARSER
  SUBDIRS += hapmixmap
endif

SUBDIRS += tools scripts doxygen

ACLOCAL_AMFLAGS = 
#use -Wno-syntax if Automake >= v1.10 to suppress "underquoted" warnings
#ACLOCAL_AMFLAGS = -Wno-syntax

EXTRA_DIST =

#create docs
docs:
	@(cd doxygen && $(MAKE) $(AM_MAKEFLAGS)) || exit 1;

