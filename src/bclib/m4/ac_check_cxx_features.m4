dnl check some of the major language features
dnl most useful for making a config.h
dnl
AC_DEFUN([AC_CHECK_CXX_FEATURES], [

OS=`uname -a`
AC_SUBST(OS)
DATE=`date`
AC_SUBST(DATE)

AH_TOP([
/******************************************************************************
 * config.h           Compiler language support flags
 *
 * This file was generated automatically when running the configure script.
 * You should rerun configure each time you switch compilers, install new
 * standard libraries, or change compiler versions.
 *
 */

])

AC_DEFINE_UNQUOTED([_compiler_name],["$CXX"],[CXX])
AC_DEFINE_UNQUOTED([_compiler_options],["$CXXFLAGS"],[CXXFLAGS])
AC_DEFINE_UNQUOTED([_os_name],["$OS"],[uname -a])
AC_DEFINE_UNQUOTED([_config_date],["$DATE"],[date])
AC_DEFINE_UNQUOTED([_platform],["$target"],[target])

AC_MSG_NOTICE([
Checking major C++ language features
])

AC_CXX_NAMESPACES
AC_CXX_EXCEPTIONS

AC_MSG_NOTICE([
Now for some of the new keywords
])
AC_CXX_TYPENAME
AC_CXX_BOOL

# support for C++0x
AX_CXX_COMPILE_STDCXX_0X

#templates
AC_CXX_TEMPLATES
#standard library

AC_MSG_NOTICE([
Checking library features
])
AC_CXX_HAVE_NUMERIC_LIMITS
AC_CXX_HAVE_STD
AC_CXX_HAVE_STL

])