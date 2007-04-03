# bcppcl, the Bayesian C++ Computation Library
# configure.ac
# Process this file with autoconf to produce a configure script.
AC_PREREQ(2.59)

AC_INIT([bcppcl],[1.0])

#put auxiliary files here
AC_CONFIG_AUX_DIR(config)

#Compute the canonical target-system type variable, target
AC_CANONICAL_TARGET

AM_INIT_AUTOMAKE([1.8 -Wall foreign])

AC_MSG_NOTICE([

Configuring $PACKAGE_STRING for $target
])

# C++ compilers
#AC_PROG_CXX
dnl automatic support for Intel, Pathscale, GNU, Portland, in that order
AC_PROG_CXX([icpc icc pathCC g++ pgCC])
AC_LANG(C++)

#preset optimisation and debug flags
AC_CXX_FLAGS_PRESET
dnl TODO: architecture-specific optimisations

# C compiler - because we are using gsl
AC_PROG_CC

# check for ranlib. remove when using libtool
AC_PROG_RANLIB

# Disable shared library building by default (use --enable-shared to enable)
#AM_DISABLE_SHARED

#check for gsl
AM_PATH_GSL

#check C++ features and write config.h
#AC_CONFIG_HEADERS(config.h)
#AC_CHECK_CXX_FEATURES

#check if script called with --enable-parallel and if so set parallel flags
AC_ENABLE_PARALLEL

# Check if the script was called with --enable-optimize, 
# --enable-debug or --enable-profile, and if so modify the compilation flags.
dnl TODO: should not allow debug and optimize together
AC_CXX_ENABLE_OPTIMIZE
AC_CXX_ENABLE_DEBUG
AC_CXX_ENABLE_PROFILE

#AC_CONFIG_SUBDIRS(utils samplers regression)

# Process the Makefile.in's
AC_CONFIG_FILES([Makefile 
	utils/Makefile 
	samplers/Makefile 
	regression/Makefile])
AC_OUTPUT

# Summarize variables
AC_ENV

#print compilation info
AC_INFO