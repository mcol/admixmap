
AC_DEFUN([AC_CXX_ENABLE_PROFILE],[

AC_MSG_CHECKING([whether to enable $CXX profiling flags])
AC_ARG_ENABLE(profile,
AS_HELP_STRING([--enable-profile],[Enable compiler profiling flags]),
[if test "$enableval" = yes; then
	AC_MSG_RESULT([yes])
  CXXFLAGS="$CXXFLAGS $CXX_PROFIL_FLAGS"
  enable_profile="yes"
fi],[
  AC_MSG_RESULT([no])
  enable_profile="no"
])
AC_ARG_VAR([CXX_PROFIL_FLAGS],[C++ compiler profiling flags])

])

