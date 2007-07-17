if PARALLEL
SUBDIRS = bclib base hapmixmap scripts doxygen
else
SUBDIRS = bclib base admixmap hapmixmap tools scripts doxygen
endif

ACLOCAL_AMFLAGS = 
#use -Wno-syntax if Automake >= v1.10 to suppress "underquoted" warnings
#ACLOCAL_AMFLAGS = -Wno-syntax

EXTRA_DIST = TODO

#create docs
docs:
	@(cd doxygen && $(MAKE) $(AM_MAKEFLAGS)) || exit 1;

