SUBDIRS		= src dist test

ACLOCAL_AMFLAGS= -I src/m4

EXTRA_DIST = INSTALL README

PKGDIR = $(top_srcdir)/packages
distcleancheck:
	@:
distuninstallcheck:
	@:
#create docs
docs:
	@(cd src/doxygen && $(MAKE) $(AM_MAKEFLAGS)) || exit 1;

#install binaries
binstall: all
	@ (cd src/admixmap && $(MAKE) $(AM_MAKEFLAGS) install) || exit 1;
	@ (cd src/hapmixmap && $(MAKE) $(AM_MAKEFLAGS) install) || exit 1;

#make admixmap binary distribution
admixmap-bindist:
	@ (cd src/base && $(MAKE) $(AM_MAKEFLAGS)) || exit 1;
	@ (cd src/admixmap && $(MAKE) $(AM_MAKEFLAGS)) || exit 1;
	rm -Rf dist/admixmap/admixmap
	mkdir -p dist/admixmap/admixmap
	cp -f src/admixmap/admixmap$(EXEEXT) dist/admixmap/admixmap
	cp -uf $(top_srcdir)/dist/admixmap/README.txt dist/admixmap/admixmap
	cp -uf $(top_srcdir)/dist/admixmap/COPYING.txt dist/admixmap/admixmap
	cp -uf $(top_srcdir)/dist/AdmixmapOutput.R dist/admixmap/admixmap
	cp -uf $(top_srcdir)/dist/admixmap/admixmap.pl dist/admixmap/admixmap
	cp -uf $(top_srcdir)/dist/doanalysis.pl dist/admixmap/admixmap
	cp -uf $(top_srcdir)/dist/admixmap/testArguments.txt dist/admixmap/admixmap
#now the tutorial
	mkdir -p dist/admixmap/admixmap/tutorial
	mkdir -p dist/admixmap/admixmap/tutorial/data
	mkdir -p dist/admixmap/admixmap/tutorial/HapFreqs
	cp -uf $(top_srcdir)/dist/admixmap/tutorial/*.pl dist/admixmap/admixmap/tutorial
	cp -uf $(top_srcdir)/dist/admixmap/tutorial/data/*.txt dist/admixmap/admixmap/tutorial/data
	cp -uf $(top_srcdir)/dist/admixmap/tutorial/HapFreqs/*.txt $(top_srcdir)/dist/admixmap/tutorial/HapFreqs/*.pl dist/admixmap/admixmap/tutorial/HapFreqs
	@echo =================================================
	@echo admixmap binary package is ready for zipping
	@echo =================================================

#make hapmixmap binary distribution
hapmixmap-bindist:
	@ (cd src/base && $(MAKE) $(AM_MAKEFLAGS)) || exit 1;
	@ (cd src/hapmixmap && $(MAKE) $(AM_MAKEFLAGS)) || exit 1;
	rm -Rf dist/hapmixmap/hapmixmap
	mkdir -p dist/hapmixmap/hapmixmap
	cp -f src/hapmixmap/hapmixmap$(EXEEXT) dist/hapmixmap/hapmixmap
	cp -uf $(top_srcdir)/dist/hapmixmap/README.txt dist/hapmixmap/hapmixmap
	cp -uf $(top_srcdir)/dist/admixmap/COPYING.txt dist/hapmixmap/hapmixmap
	cp -uf $(top_srcdir)/dist/AdmixmapOutput.R dist/hapmixmap/hapmixmap
#	cp -uf $(top_srcdir)/test/hapmixmap/hapmixmap.pl dist/hapmixmap/hapmixmap
	cp -uf $(top_srcdir)/dist/doanalysis.pl dist/hapmixmap/hapmixmap
	cp -uf $(top_srcdir)/dist/hapmixmap/options.txt dist/hapmixmap/hapmixmap
	cp -uf $(top_srcdir)/dist/hapmixmap/tutorial.pl dist/hapmixmap/hapmixmap
	cp -uf $(top_srcdir)/dist/hapmixmap/tutorial.html dist/hapmixmap/hapmixmap
	cp -uf $(top_srcdir)/dist/hapmixmap/*.conf dist/hapmixmap/hapmixmap
	mkdir dist/hapmixmap/hapmixmap/data
	cp -uf $(top_srcdir)/dist/hapmixmap/data/*.txt dist/hapmixmap/hapmixmap/data
	@echo ==================================================
	@echo hapmixmap binary package is ready for zipping
	@echo ==================================================

#make admixmap (only) source distribution
admixmap-srcdist: dist
#unpack dist
	tar -xzf @PACKAGE_NAME@-@PACKAGE_VERSION@.tar.gz
#copy sources
	cp -Rf @PACKAGE_NAME@-@PACKAGE_VERSION@/src/admixmap $(PKGDIR)/admixmap
	cp -Rf @PACKAGE_NAME@-@PACKAGE_VERSION@/src/base $(PKGDIR)/admixmap
	cp -Rf @PACKAGE_NAME@-@PACKAGE_VERSION@/src/m4 $(PKGDIR)/admixmap
#copy other files
	cp -Rf @PACKAGE_NAME@-@PACKAGE_VERSION@/dist/admixmap/*.txt $(PKGDIR)/admixmap
	cp -uf @PACKAGE_NAME@-@PACKAGE_VERSION@/INSTALL $(PKGDIR)/admixmap
	cp -uf @PACKAGE_NAME@-@PACKAGE_VERSION@/dist/AdmixmapOutput.R $(PKGDIR)/admixmap
	cp -uf @PACKAGE_NAME@-@PACKAGE_VERSION@/dist/admixmap/admixmap.pl $(PKGDIR)/admixmap
	cp -uf @PACKAGE_NAME@-@PACKAGE_VERSION@/dist/doanalysis.pl $(PKGDIR)/admixmap
	cp -uf @PACKAGE_NAME@-@PACKAGE_VERSION@/dist/admixmap/testArguments.txt $(PKGDIR)/admixmap
	cp -rf @PACKAGE_NAME@-@PACKAGE_VERSION@/dist/admixmap/tutorial $(PKGDIR)/admixmap
#delete temporaries
	rm -Rf @PACKAGE_NAME@-@PACKAGE_VERSION@
	rm -f @PACKAGE_NAME@-@PACKAGE_VERSION@.tar.gz
#build package
	cd $(PKGDIR)/admixmap \
&& autoreconf -fi \
&& ./configure \
&& make dist
	mv $(PKGDIR)/admixmap/admixmap-*.tar.gz .
	@echo =================================================
	@echo admixmap source package is ready for distribution
	@echo =================================================

#make hapmixmap (only) source distribution
hapmixmap-srcdist: dist
#unpack dist
	tar -xzf @PACKAGE_NAME@-@PACKAGE_VERSION@.tar.gz
#copy sources
	cp -Rf @PACKAGE_NAME@-@PACKAGE_VERSION@/src/hapmixmap $(PKGDIR)/hapmixmap
	cp -Rf @PACKAGE_NAME@-@PACKAGE_VERSION@/src/base $(PKGDIR)/hapmixmap
	cp -Rf @PACKAGE_NAME@-@PACKAGE_VERSION@/src/m4 $(PKGDIR)/hapmixmap
#formatter files
	cp -Rf @PACKAGE_NAME@-@PACKAGE_VERSION@/src/tools/hapmixmapFormatter $(PKGDIR)/hapmixmap
#copy other files
	cp -uf @PACKAGE_NAME@-@PACKAGE_VERSION@/dist/hapmixmap/README.txt $(PKGDIR)/hapmixmap
	cp -uf @PACKAGE_NAME@-@PACKAGE_VERSION@/dist/admixmap/COPYING.txt $(PKGDIR)/hapmixmap
	cp -uf @PACKAGE_NAME@-@PACKAGE_VERSION@/INSTALL $(PKGDIR)/hapmixmap
	cp -uf @PACKAGE_NAME@-@PACKAGE_VERSION@/dist/AdmixmapOutput.R $(PKGDIR)/hapmixmap
#	cp -uf @PACKAGE_NAME@-@PACKAGE_VERSION@/dist/hapmixmap/hapmixmap.pl $(PKGDIR)/hapmixmap
	cp -uf @PACKAGE_NAME@-@PACKAGE_VERSION@/dist/doanalysis.pl $(PKGDIR)/hapmixmap
	cp -uf @PACKAGE_NAME@-@PACKAGE_VERSION@/dist/hapmixmap/options.txt $(PKGDIR)/hapmixmap
	cp -uf @PACKAGE_NAME@-@PACKAGE_VERSION@/dist/hapmixmap/tutorial.pl $(PKGDIR)/hapmixmap
	cp -uf @PACKAGE_NAME@-@PACKAGE_VERSION@/dist/hapmixmap/tutorial.html $(PKGDIR)/hapmixmap
	cp -Rf @PACKAGE_NAME@-@PACKAGE_VERSION@/dist/hapmixmap/*.conf $(PKGDIR)/hapmixmap
	cp -Rf @PACKAGE_NAME@-@PACKAGE_VERSION@/dist/hapmixmap/data $(PKGDIR)/hapmixmap
#delete temporaries
	rm -Rf @PACKAGE_NAME@-@PACKAGE_VERSION@
	rm -f @PACKAGE_NAME@-@PACKAGE_VERSION@.tar.gz
#build package
	cd $(PKGDIR)/hapmixmap \
&& autoreconf -fi \
&& ./configure \
&& make dist
	mv $(PKGDIR)/hapmixmap/hapmixmap-*.tar.gz .
	@echo ==================================================
	@echo hapmixmap source package is ready for distribution
	@echo ==================================================

#make bcppcl (only) source distribution
bcppcl-srcdist:
	@(cd $(top_srcdir)/src/bcppcl && autoreconf -fi)
	@($(top_srcdir)/src/bcppcl/configure && \
make dist) || exit 1;
	@echo ===============================================
	@echo bcppcl source package is ready for distribution
	@echo ===============================================
