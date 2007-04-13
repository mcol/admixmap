AC_DEFUN([AC_SPRNG], [

AC_ARG_WITH(sprng-prefix,
[--with-sprng-prefix=PFX   Prefix where SPRNG is installed (optional)],
            sprng_prefix="$withval", sprng_prefix="")

if test x = x"$SPRNGLIBS"; then
dnl TODO: check for sprng lib. The following does not work.
dnl  AC_CHECK_FUNC(sprng, [SPRNGLIBS=" "])
dnl  AC_CHECK_LIB(sprng, init_sprng, [SPRNGLIBS="-lsprng"])
  AC_CHECK_HEADER([sprng.h],
[SPRNGLIBS="-lsprng"], [AC_MSG_ERROR([Cannot find sprng.h])])
fi

AC_SUBST(SPRNGLIBS)

# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x = x"$SPRNGLIBS"; then
        $2
        :
else
        ifelse([$1],,[AC_DEFINE(HAVE_SPRNG,1,[Define if you have the SPRNG library.])],[$1])
        :
fi
])dnl AC_SPRNG

