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

#make admixmap binary distribution
admixmap-bindist:
	@ (cd src/admixmap && $(MAKE) $(AM_MAKEFLAGS)) || exit 1;
	rm -Rf dist/admixmap/admixmap
	mkdir -p dist/admixmap/admixmap
	cp -f src/admixmap/admixmap$(EXEEXT) dist/admixmap/admixmap
	cp -uf dist/admixmap/README.txt dist/admixmap/admixmap
	cp -uf dist/admixmap/COPYING.txt dist/admixmap/admixmap
	cp -uf test/admixmap/AdmixmapOutput.R dist/admixmap/admixmap
	cp -uf test/admixmap/admixmap.pl dist/admixmap/admixmap
	cp -uf dist/doanalysis.pl dist/admixmap/admixmap
	cp -uf test/admixmap/testArguments.txt dist/admixmap/admixmap
	cp -rf test/admixmap/tutorial dist/admixmap/admixmap
	@( cd dist/admixmap \
&& rm -Rf admixmap/tutorial/.svn \
&& rm -f admixmap-bin.tar.gz \
&& tar -czf admixmap-bin.tar.gz admixmap \
&& rm -Rf admixmap) || exit 1;

#make hapmixmap binary distribution
hapmixmap-bindist:
	@ (cd src/hapmixmap && $(MAKE) $(AM_MAKEFLAGS)) || exit 1;
	rm -Rf dist/hapmixmap/hapmixmap
	mkdir -p dist/hapmixmap/hapmixmap
	cp -f src/hapmixmap/hapmixmap$(EXEEXT) dist/hapmixmap/hapmixmap
	cp -uf dist/hapmixmap/README.txt dist/hapmixmap/hapmixmap
	cp -uf dist/admixmap/COPYING.txt dist/hapmixmap/hapmixmap
	cp -uf test/admixmap/AdmixmapOutput.R dist/hapmixmap/hapmixmap
#	cp -uf test/hapmixmap/hapmixmap.pl dist/hapmixmap/hapmixmap
	cp -uf dist/doanalysis.pl dist/hapmixmap/hapmixmap
	cp -uf dist/hapmixmap/options.txt dist/hapmixmap/hapmixmap
	cp -uf dist/hapmixmap/tutorial.pl dist/hapmixmap/hapmixmap
	cp -uf dist/hapmixmap/tutorial.html dist/hapmixmap/hapmixmap
#	cp -uf dist/hapmixmap/*.conf dist/hapmixmap/hapmixmap
	mkdir dist/hapmixmap/hapmixmap/data
	cp -Rf dist/hapmixmap/data/*.txt dist/hapmixmap/hapmixmap/data
	@( cd dist/hapmixmap \
&& rm -f hapmixmap-bin.tar.gz \
&& tar -czf hapmixmap-bin.tar.gz hapmixmap \
&& rm -Rf hapmixmap) || exit 1;

#make admixmap (only) source distribution
admixmap-srcdist: dist
#unpack dist
	tar -xzf @PACKAGE_NAME@-@PACKAGE_VERSION@.tar.gz
#copy sources
	cp -Rf @PACKAGE_NAME@-@PACKAGE_VERSION@/src/admixmap packages/admixmap
	cp -Rf @PACKAGE_NAME@-@PACKAGE_VERSION@/src/base packages/admixmap
	cp -Rf @PACKAGE_NAME@-@PACKAGE_VERSION@/src/m4 packages/admixmap
#copy other files
	cp -Rf @PACKAGE_NAME@-@PACKAGE_VERSION@/dist/admixmap/*.txt packages/admixmap
	cp -uf @PACKAGE_NAME@-@PACKAGE_VERSION@/INSTALL packages/admixmap
	cp -uf @PACKAGE_NAME@-@PACKAGE_VERSION@/test/admixmap/AdmixmapOutput.R packages/admixmap
	cp -uf @PACKAGE_NAME@-@PACKAGE_VERSION@/test/admixmap/admixmap.pl packages/admixmap
	cp -uf @PACKAGE_NAME@-@PACKAGE_VERSION@/dist/doanalysis.pl packages/admixmap
	cp -uf @PACKAGE_NAME@-@PACKAGE_VERSION@/test/admixmap/testArguments.txt packages/admixmap
	cp -rf @PACKAGE_NAME@-@PACKAGE_VERSION@/test/admixmap/tutorial packages/admixmap
#delete temporaries
	rm -Rf @PACKAGE_NAME@-@PACKAGE_VERSION@
	rm -f @PACKAGE_NAME@-@PACKAGE_VERSION@.tar.gz
#build package
	cd packages/admixmap \
&& autoreconf -fi \
&& ./configure \
&& make dist
	mv packages/admixmap/admixmap-*.tar.gz .
	@echo =================================================
	@echo admixmap source package is ready for distribution
	@echo =================================================
## TODO: check binary path

#make hapmixmap (only) source distribution
hapmixmap-srcdist: dist
#unpack dist
	tar -xzf @PACKAGE_NAME@-@PACKAGE_VERSION@.tar.gz
#copy sources
	cp -Rf @PACKAGE_NAME@-@PACKAGE_VERSION@/src/hapmixmap packages/hapmixmap
	cp -Rf @PACKAGE_NAME@-@PACKAGE_VERSION@/src/base packages/hapmixmap
	cp -Rf @PACKAGE_NAME@-@PACKAGE_VERSION@/src/m4 packages/hapmixmap
#copy other files
	cp -uf @PACKAGE_NAME@-@PACKAGE_VERSION@/dist/hapmixmap/README.txt packages/hapmixmap
	cp -uf @PACKAGE_NAME@-@PACKAGE_VERSION@/dist/admixmap/COPYING.txt packages/hapmixmap
	cp -uf @PACKAGE_NAME@-@PACKAGE_VERSION@/INSTALL packages/hapmixmap
	cp -uf @PACKAGE_NAME@-@PACKAGE_VERSION@/test/admixmap/AdmixmapOutput.R packages/hapmixmap
#	cp -uf @PACKAGE_NAME@-@PACKAGE_VERSION@/test/hapmixmap/hapmixmap.pl packageshapmixmap
	cp -uf @PACKAGE_NAME@-@PACKAGE_VERSION@/dist/doanalysis.pl packages/hapmixmap
	cp -uf @PACKAGE_NAME@-@PACKAGE_VERSION@/dist/hapmixmap/options.txt packages/hapmixmap
	cp -uf @PACKAGE_NAME@-@PACKAGE_VERSION@/dist/hapmixmap/tutorial.pl packages/hapmixmap
	cp -uf @PACKAGE_NAME@-@PACKAGE_VERSION@/dist/hapmixmap/tutorial.html packages/hapmixmap
#	cp -uf @PACKAGE_NAME@-@PACKAGE_VERSION@/dist/hapmixmap/*.conf packages/hapmixmap
	cp -Rf @PACKAGE_NAME@-@PACKAGE_VERSION@/dist/hapmixmap/data packages/hapmixmap
#delete temporaries
	rm -Rf @PACKAGE_NAME@-@PACKAGE_VERSION@
	rm -f @PACKAGE_NAME@-@PACKAGE_VERSION@.tar.gz
#build package
	cd packages/hapmixmap \
&& autoreconf -fi \
&& ./configure \
&& make dist
	mv packages/hapmixmap/hapmixmap-*.tar.gz .
	@echo ==================================================
	@echo hapmixmap source package is ready for distribution
	@echo ==================================================

#make bcppcl (only) source distribution
bcppcl-srcdist:
	@(cd src/bcppcl && \
autoreconf -fi && \
./configure && \
make dist) || exit 1;
	mv src/bcppcl/bcppcl-*.tar.gz .
	@echo ===============================================
	@echo bcppcl source package is ready for distribution
	@echo ===============================================
