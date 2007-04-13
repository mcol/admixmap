AC_PREREQ(2.59)

AC_INIT([genepi], [x])

#put auxiliary files here
AC_CONFIG_AUX_DIR(config)

AC_ARG_ENABLE(bcppcl,
AS_HELP_STRING([--enable-bcppcl],[Enable building of bcppcl (yes)]),
[bcppcl_enabled=yes],[bcppcl_enabled=yes])

if test "$bcppcl_enabled" = yes; then
  AC_CONFIG_SUBDIRS([bcppcl])
fi

AM_INIT_AUTOMAKE([1.8 -Wall foreign])
AC_CONFIG_FILES([Makefile 
])

AC_OUTPUT