lib_LTLIBRARIES	= libbayes.la
libbayes_la_SOURCES	= bayes.cc
libbayes_la_LIBADD	= samplers/libsampler.la			\
			  utils/libutils.la				\
			  regression/libbayesreg.la
SUBDIRS			= samplers utils regression
noinst_HEADERS		= common.h bcppcl.h

if PARALLEL
lib_LTLIBRARIES			+= libbayes-para.la
libbayes_para_la_SOURCES	= $(libbayes_la_SOURCES)
libbayes_para_la_LIBADD		= samplers/libsampler-para.la		\
				  utils/libutils-para.la		\
				  regression/libbayesreg-para.la
libbayes_para_CPPFLAGS		= -DPARALLEL
endif

# libbayes_a_LDFLAGS	= -version-info 0:0:0
