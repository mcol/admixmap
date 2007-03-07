# noinst_PROGRAMS			= rs-lookup
AM_YFLAGS = -d
noinst_PROGRAMS		= 

#SUBDIRS = DataFormatter
 
# Programs and scripts to be installed in the $prefix/bin directory
bin_PROGRAMS		= 

bin_SCRIPTS		=

noinst_HEADERS		= Legend.h \
			  UserData.h \
			  UserDataIndividual.h \
			  UserDataLocus.h \
			  DataFormatter.h \
			  HapmixmapDataFormatter.h \
			  ge-snp-lookup.h \
			  itoa.h

man_MANS		= 

BUILT_SOURCES =

if CONVERTER
bin_PROGRAMS		+= ge-snp-lookup ge-user2hmx
bin_SCRIPTS		+= ge-build-snp-list ge-fetch-legend
man_MANS		+= ge-user2hmx.1
ge_snp_lookup_SOURCES	= ge-snp-lookup.cc \
			  Legend.cc

ge_user2hmx_SOURCES	= \
				  UserData.cc \
				  UserDataIndividual.cc \
				  DataFormatter.cc \
				  HapmixmapDataFormatter.cc \
				  UserDataLocus.cc \
				  Legend.cc \
				  userdata-lexer.ll \
				  userdata-parser.yy

BUILT_SOURCES 		+= userdata-parser.h
endif

#
# Information about how to build the programs
#

EXTRA_DIST		= $(man_MANS) $(bin_SCRIPTS)
