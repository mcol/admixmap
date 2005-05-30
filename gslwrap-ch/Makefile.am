## Makefile.am -- Process this file with automake to produce Makefile.in
## Copyright (C) 2000 Gary V. Vaughan
##
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2, or (at your option)
## any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; if not, write to the Free Software
## Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

libtool: @LIBTOOL_DEPS@
	cd $(srcdir) && \
		$(SHELL) ./config.status --recheck

LIBS = -lm -lgsl -lgslcblas

MAINTAINERCLEANFILES	= Makefile.in

INCLUDES		= -I$(top_builddir)
pkginclude_HEADERS 	= matrix.h vector.h matrix_d.h vector_d.h \
											matrix_i.h vector_i.h \
											rand.h 

lib_LTLIBRARIES 		= libmat.la

libmat_la_SOURCES	= matrix.cc vector.cc matrix_d.cc vector_d.cc \
											matrix_i.cc vector_i.cc rand_serial.cc

libmat_la_LDFLAGS = -static
