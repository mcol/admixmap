SUBDIRS		= src dist test

ACLOCAL_AMFLAGS= -I src/m4

EXTRA_DIST = INSTALL README

distcleancheck:
	@:
distuninstallcheck:
	@:
#create docs
docs:
	@(cd src/doxygen && $(MAKE) $(AM_MAKEFLAGS)) || exit 1;
