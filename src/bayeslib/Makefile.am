noinst_LTLIBRARIES	= libbayes.la
libbayes_la_SOURCES	= bayes.cc
# libbayes_a_LDFLAGS	= -version-info 0:0:0
libbayes_la_LIBADD	= samplers/libsampler.la \
			  utils/libutils.la \
			  regression/libbayesreg.la
SUBDIRS			= samplers utils regression
noinst_HEADERS		= common.h
