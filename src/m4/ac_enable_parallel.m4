
AC_DEFUN([AC_ENABLE_PARALLEL],
[
AC_MSG_CHECKING([whether to build the parallel version])
AC_ARG_ENABLE(parallel, 
AS_HELP_STRING([--enable-parallel],[Build the parallel (distributed-memory/OpenMPI) version.
Contrast with --enable-openmp which adds shared-memory parallelism.]),
[if test "$enableval" = yes; then
	AC_MSG_RESULT([yes])
	enable_parallel="yes"
        #AC_REQUIRE([ACX_MPI])
	PARALLEL_ENABLE="yes"
fi],[
  AC_MSG_RESULT([no])
  enable_parallel="no"
  PARALLEL_ENABLE="no"
])

AM_CONDITIONAL(PARALLEL, test x$enable_parallel = xyes)

if test "x$enable_parallel" = "xyes";
  then AC_SPRNG AC_CXX_MPI
       PARALLEL_LIBS="$MPILIBS -lgmp $SPRNGLIBS"
fi
AC_SUBST([PARALLEL_ENABLE])
AC_SUBST([PARALLEL_LIBS])


])
