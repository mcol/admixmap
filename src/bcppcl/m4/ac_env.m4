
AC_DEFUN([AC_ENV],[
AC_MSG_NOTICE([
------------------------------------------------------------------------
Configuration:

C Compiler         = $CC
C++ Compiler       = $CXX
Compiler flags     = $CXXFLAGS
dnl CXX_OPTIMIZE_FLAGS = $CXX_OPTIMIZE_FLAGS
dnl CXX_DEBUG_FLAGS    = $CXX_DEBUG_FLAGS
Optimize           = $enable_optimize
debug              = $enable_debug
profile            = $enable_profile
Parallel version   = $enable_parallel
Host System        = $host
Install path:      = $prefix
------------------------------------------------------------------------
])])

