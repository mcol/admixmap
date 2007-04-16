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

admixmap-bindist:
	@ (cd src/admixmap && $(MAKE) $(AM_MAKEFLAGS)) || exit 1;
	rm -Rf dist/admixmap/admixmap
	mkdir -p dist/admixmap/admixmap
	cp -f src/admixmap/admixmap$(EXEEXT) dist/admixmap/admixmap
	cp -uf dist/admixmap/README.txt dist/admixmap/admixmap
	cp -uf dist/admixmap/COPYING.txt dist/admixmap/admixmap
	cp -uf test/admixmap/AdmixmapOutput.R dist/admixmap/admixmap
	cp -uf test/admixmap/admixmap.pl dist/admixmap/admixmap
	cp -uf test/admixmap/testArguments.txt dist/admixmap/admixmap
	cp -rf test/admixmap/tutorial dist/admixmap/admixmap
	@( cd dist/admixmap \
&& rm -Rf admixmap/tutorial/.svn \
&& rm -f admixmap-bin.tar.gz \
&& tar -czf admixmap-bin.tar.gz admixmap \
&& rm -Rf admixmap) || exit 1;

hapmixmap-bindist:
	@ (cd src/hapmixmap && $(MAKE) $(AM_MAKEFLAGS)) || exit 1;
	rm -Rf dist/hapmixmap/hapmixmap
	mkdir -p dist/hapmixmap/hapmixmap
	cp -f src/hapmixmap/hapmixmap$(EXEEXT) dist/hapmixmap/hapmixmap
	cp -uf dist/hapmixmap/README.txt dist/hapmixmap/hapmixmap
	cp -uf dist/admixmap/COPYING.txt dist/hapmixmap/hapmixmap
	cp -uf test/admixmap/AdmixmapOutput.R dist/hapmixmap/hapmixmap
#	cp -uf test/hapmixmap/hapmixmap.pl dist/hapmixmap/hapmixmap
	cp -uf dist/hapmixmap/options.txt dist/hapmixmap/hapmixmap
	cp -uf dist/hapmixmap/tutorial.pl dist/hapmixmap/hapmixmap
	cp -uf dist/hapmixmap/tutorial.html dist/hapmixmap/hapmixmap
	cp -uf dist/hapmixmap/*.conf dist/hapmixmap/hapmixmap
	mkdir dist/hapmixmap/hapmixmap/data
	cp -Rf dist/hapmixmap/data/*.txt dist/hapmixmap/hapmixmap/data
	@( cd dist/hapmixmap \
&& rm -f hapmixmap-bin.tar.gz \
&& tar -czf hapmixmap-bin.tar.gz hapmixmap \
&& rm -Rf hapmixmap) || exit 1;

