AC_PREREQ(2.59)

AC_INIT([bclib],[1.2])

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
AC_PROG_CXX([icpc pathCC g++ pgCC])
AC_PROG_CC([icc pathcc gcc pgcc])
AC_LANG(C++)

#preset optimisation and debug flags
AC_CXX_FLAGS_PRESET
dnl TODO: architecture-specific optimisations

# C compiler - because we are using gsl
AC_PROG_CC

# check for ranlib. not stricly necessary but just in case 
#AC_PROG_RANLIB

#check C++ features
#AC_CONFIG_HEADERS(config.h)
#AC_CHECK_CXX_FEATURES

#check for std headers
#AC_CHECK_HEADERS(vector,,AC_MSG_WARN(STL classes missing ?))
#AC_CHECK_HEADERS(string,,AC_MSG_WARN(STL classes missing ?))

#check for gsl
AM_PATH_GSL

# check if script called with --enable-parallel and if so,
# check for mpi compiler wrapper, MPI library and MPI header,
# then set parallel flags and change CXX to the MPI equivalent
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

#AC_CONFIG_SUBDIRS(utils samplers regression)

# Process the Makefile.in's
AC_CONFIG_FILES([Makefile 
	utils/Makefile 
	samplers/Makefile 
	regression/Makefile
        bclib/Makefile
        examples/Makefile
	doxygen/Makefile
	doxygen/Doxyfile
	])
AC_OUTPUT

# Summarize variables
AC_ENV

#print compilation info
AC_INFO

