
AC_DEFUN([AC_ENV],[
AC_MSG_NOTICE([
------------------------------------------------------------------------
Configuration:

C Compiler              = $CC
C++ Compiler            = $CXX
Compiler flags          = $WFLAGS $CXXFLAGS
dnl CXX_OPTIMIZE_FLAGS  = $CXX_OPTIMIZE_FLAGS
dnl CXX_DEBUG_FLAGS     = $CXX_DEBUG_FLAGS
Optimize                = $enable_optimize
debug                   = $enable_debug
Range-checking          = ${enable_range_check:-no}
OpenMP parallelization  = ${enable_openmp:-no}
profile                 = $enable_profile
Platform                = $target
Install path:           = $prefix
------------------------------------------------------------------------
])

])

