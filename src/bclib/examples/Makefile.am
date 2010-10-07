noinst_PROGRAMS = OptionReaderDemo RObjectDemo

AM_CPPFLAGS = -I$(srcdir)/..
AM_LDFLAGS = -L$(top_builddir)/$(subdir)/..
LDADD = -lbcl
AM_CXXFLAGS=$(WFLAGS) -std=c++0x
SUFFIXES=.cpp .o
OptionReaderDemo_SOURCES = OptionReaderDemo.cpp
RObjectDemo_SOURCES = RObjectDemo.cpp