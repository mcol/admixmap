AC_PREREQ(2.59)

AC_INIT([genepi],[src])

#put auxiliary files here
AC_CONFIG_AUX_DIR(config)

#Compute the canonical target-system type variable, target
AC_CANONICAL_TARGET

dnl -Wno-override is required because some default makefile
dnl targets are overridden to make 'make distcheck' work
AM_INIT_AUTOMAKE([1.8 -Wall -Wno-override foreign])

AC_MSG_NOTICE([

Configuring $PACKAGE_STRING for $target
])

# C++ compilers
#AC_PROG_CXX
AC_PROG_CXX([icpc icc pathCC g++ pgCC])
AC_LANG(C++)

#preset optimisation and debug flags
AC_CXX_FLAGS_PRESET
dnl TODO: architecture-specific optimisations

# C compiler - because we are using gsl
AC_PROG_CC

#check C++ features
AC_CONFIG_HEADERS(config.h)
AC_CHECK_CXX_FEATURES

#check for std headers
#AC_CHECK_HEADERS(vector,,AC_MSG_WARN(STL classes missing ?))
#AC_CHECK_HEADERS(string,,AC_MSG_WARN(STL classes missing ?))

#check for gsl
AM_PATH_GSL

#check if script called with --enable-parallel and if so set parallel flags
AC_ENABLE_PARALLEL

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
#AC_PROG_RANLIB

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
	   src/scripts/Makefile
           src/scripts/genepi-config
	   src/tools/Makefile
           src/tools/hapmixmapFormatter/Makefile
           src/tools/hapmix2impute/Makefile
           src/tools/ancestrymapConversion/Makefile
	   src/doxygen/Makefile
	   src/doxygen/Doxyfile
	   test/Makefile
	   test/admixmap/Makefile
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

