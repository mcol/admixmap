noinst_LTLIBRARIES         	= libbayesreg.la

libbayesreg_la_SOURCES		= Regression.cc \
				  LinearRegression.cc \
				  LogisticRegression.cc \
				  CoxRegression.cc

noinst_HEADERS			= Regression.h \
				  LinearRegression.h \
				  LogisticRegression.h \
				  CoxRegression.h

INCLUDES			= -I$(top_srcdir)/src/bayeslib

if PARALLEL
noinst_LTLIBRARIES		+= libbayesreg-para.la
libbayesreg_para_la_SOURCES	= $(libbayesreg_la_SOURCES)
libbayesreg_para_la_CPPFLAGS	= -DPARALLEL
endif

# libbayesreg_la_LDFLAGS		= -version-info 0:0:0
