AC_PREREQ(2.59)

AC_INIT([genepi],[src])

#put auxiliary files here
AC_CONFIG_AUX_DIR(config)

## DDF: added this to eliminate warning from 'autoreconf' regarding libtool:
#put additional local autoconf macros here
AC_CONFIG_MACRO_DIR([src/m4])

#Compute the canonical target-system type variable, target
AC_CANONICAL_TARGET

dnl -Wno-override is required because some default makefile
dnl targets are overridden to make 'make distcheck' work
AM_INIT_AUTOMAKE([1.8 -Wall -Wno-override foreign])

AC_MSG_NOTICE([

Configuring $PACKAGE_STRING for $target
])


# --enable-old-parser which is used to set or clear USE_GENOTYPE_PARSER
AC_ARG_ENABLE(old-parser,
	[  --enable-old-parser     don't use the new GFileLexer/SimpleLocusParser/GenotypeParser
			  classes to parse input files (required for pedfiles)])
AM_CONDITIONAL(USE_GENOTYPE_PARSER,test "$enable_old_parser" != yes)
if test "$enable_old_parser" = yes; then
    dnl AC_SUBST(USE_GENOTYPE_PARSER,0)
    AC_DEFINE([USE_GENOTYPE_PARSER],[0],[switch which input-file parser to use])
else
    dnl AC_SUBST(USE_GENOTYPE_PARSER,1)
    AC_DEFINE([USE_GENOTYPE_PARSER],[1],[switch which input-file parser to use])
fi


# --enable-range-check which will set or clear AGGRESSIVE_RANGE_CHECK
# Consider enabling by default for --enable-debug:
#if test "$disable_range_check" = no -o "$enable_debug" = yes ; then
AC_ARG_ENABLE(range-check,
	[  --enable-range-check   turn on internal range-checking of array indices
			  (very useful for debugging)])
if test "$enable_range_check" = yes ; then
    AC_DEFINE([AGGRESSIVE_RANGE_CHECK],[1],[switch on/off array index range checking])
else
    AC_DEFINE([AGGRESSIVE_RANGE_CHECK],[0],[switch on/off array index range checking])
fi


# Produce a less verbose output from the makefiles with automake 1.11.
# To produce the usual output, run 'make V=1'.
m4_ifdef([AM_SILENT_RULES], [AM_SILENT_RULES([yes])])


# C++ compilers
#AC_PROG_CXX
AC_PROG_CXX([icpc icc pathCC g++ pgCC openCC])
AC_LANG(C++)

#preset optimisation and debug flags
AC_CXX_FLAGS_PRESET
dnl TODO: architecture-specific optimisations

#check C++ features
AC_CONFIG_HEADERS(config.h)
AC_CHECK_CXX_FEATURES

#check for std headers
#AC_CHECK_HEADERS(vector,,AC_MSG_WARN(STL classes missing ?))
#AC_CHECK_HEADERS(string,,AC_MSG_WARN(STL classes missing ?))

#check for gsl
AM_PATH_GSL

#check for boost (optional)
AC_CHECK_HEADERS([boost.h])

#check if script called with --enable-parallel and if so set parallel flags
AC_ENABLE_PARALLEL

# --enable-openmp
AC_ARG_ENABLE(openmp,
AS_HELP_STRING([--enable-openmp],[Use OpenMP to parallelize computations across
multiple shared-memory processors.  Contrast with --enable-parallel which uses
OpenMPI to parallelize across multiple distributed-memory nodes.]))
AM_CONDITIONAL(USE_OPENMP,test "$enable_openmp" = yes)


# Check if the script was called with --enable-optimize,
# --enable-debug or --enable-profile, and if so modify the compilation flags.
dnl TODO: should not allow debug and optimize together
AC_CXX_ENABLE_OPTIMIZE
AC_CXX_ENABLE_DEBUG
AC_CXX_ENABLE_PROFILE

# Doxygen documentation
AC_PROG_DOXYGEN

# Disable shared library building by default (use --enable-shared to enable)
AM_DISABLE_SHARED
#invoke libtool - must be *after* enable checks
AC_PROG_LIBTOOL
LT_INIT

## write definitions of compiler, flags etc in config.h
AC_DEFINE_UNQUOTED([_compiler_name],["$CXX"],[C++ compiler name])
AC_DEFINE_UNQUOTED([_compiler_flags],["$CXXFLAGS"],[Compiler flags])
AC_DEFINE_UNQUOTED([_os_name],["$OS"],[uname -a])
AC_DEFINE_UNQUOTED([_config_date],["$DATE"],[date])
AC_DEFINE_UNQUOTED([_host],["$host"],[Built on])
AC_DEFINE_UNQUOTED([_platform],["$target"],[target platform])

# Process the Makefile.in's
AC_CONFIG_FILES([Makefile
	   src/Makefile
	   src/admixmap/Makefile
	   src/hapmixmap/Makefile
	   src/base/Makefile
	   src/bclib/Makefile
	   src/bclib/bclib/Makefile
	   src/bclib/regression/Makefile
	   src/bclib/samplers/Makefile
	   src/bclib/utils/Makefile
	   src/bclib/examples/Makefile
	   src/bclib/doxygen/Makefile
	   src/bclib/doxygen/Doxyfile
	   src/scripts/Makefile
	   src/scripts/genepi-config
	   src/tools/Makefile
	   src/tools/hapmap/Makefile
	   src/tools/FPHD/Makefile
	   src/tools/hapmix2impute/Makefile
	   src/tools/ancestrymapConversion/Makefile
	   src/tools/bir/Makefile
	   src/doxygen/Makefile
	   src/doxygen/Doxyfile
	   test/Makefile
	   test/admixmap/Makefile
	   test/admixmap/pedigree-tester/Makefile
	   test/hapmixmap/Makefile
	   dist/Makefile
	   dist/admixmap/Makefile
	   dist/admixmap/tutorial/Makefile
	   dist/admixmap/tutorial/data/Makefile
	   dist/admixmap/tutorial/HapFreqs/Makefile
	   dist/hapmixmap/Makefile
	   dist/hapmixmap/data/Makefile
	])
AC_OUTPUT

# Summarize variables
AC_ENV

#print compilation info
AC_INFO
