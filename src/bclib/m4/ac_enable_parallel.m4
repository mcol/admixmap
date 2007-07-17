
AC_DEFUN([AC_ENABLE_PARALLEL],
[
AC_MSG_CHECKING([whether to build the parallel version])
AC_ARG_ENABLE(parallel, 
AS_HELP_STRING([--enable-parallel],[build the parallel version]),
[if test "$enableval" = yes; then
	AC_MSG_RESULT([yes])
	enable_parallel="yes"
	#AC_CHECK_HEADERS([mpi++.h],,AC_MSG_ERROR(cannot find MPI headers))
        #AC_REQUIRE([ACX_MPI])
	PARALLEL_ENABLE="yes"
fi],[
  AC_MSG_RESULT([no])
  enable_parallel="no"
  PARALLEL_ENABLE="no"
])

AM_CONDITIONAL(PARALLEL, test x$enable_parallel = xyes)

if test "x$enable_parallel" = "xyes";
  then AC_CXX_MPI
       PARALLEL_LIBS="$MPILIBS -lgmp -lsprng"
fi
AC_SUBST([PARALLEL_ENABLE])
AC_SUBST([PARALLEL_LIBS])


])
