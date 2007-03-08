noinst_LTLIBRARIES		= libinterface.la

libinterface_la_SOURCES		= IInputData.cc \
				  IInputData.h \
				  IOptions.cc \
				  IOptions.h \
				  IHapMixFreqs.cc \
				  IHapMixFreqs.h \
				  IGenome.cc \
				  IGenome.h \
				  IFreqArray.cc \
				  IFreqArray.h \
				  ICompositeLocus.cc \
				  ICompositeLocus.h \
				  IChromosome.cc \
				  IChromosome.h

INCLUDES			= -I$(top_srcdir)/src/bayeslib

if PARALLEL
noinst_LTLIBRARIES		+= libinterface-para.la
libinterface_para_la_SOURCES	= $(libinterface_la_SOURCES)
libinterface_para_la_CPPFLAGS	= -DPARALLEL
endif
