
dnl Compiler flags preset
AC_DEFUN([AC_CXX_FLAGS_PRESET],[

dnl Declare variables which we want substituted in the Makefile.in's

dnl AC_SUBST(CXXFLAGS)
AC_SUBST(CXX_OPTIMIZE_FLAGS)
AC_SUBST(CXX_DEBUG_FLAGS)
AC_SUBST(CXX_PROFIL_FLAGS)
dnl AC_SUBST(CXX_LIBS)
AC_SUBST(AR)
AC_SUBST(AR_FLAGS)
dnl AC_SUBST(LDFLAGS)
AC_SUBST(RANLIB)

dnl Set default values
AR=ar
AR_FLAGS="-cru"

AC_ARG_VAR([EXTRA_FLAGS], [Extra C++ compiler flags])

AC_MSG_CHECKING([whether using $CXX preset flags])
AC_ARG_ENABLE(cxx-flags-preset,
AS_HELP_STRING([--enable-cxx-flags-preset],
[Enable C++ compiler flags preset @<:@default yes@:>@]),[],[enableval='yes'])

if test "$enableval" = yes ; then
        CXXFLAGS=$EXTRA_FLAGS
	ac_cxx_flags_preset=yes

	case "$CXX" in
	*icpc*|*icc*) dnl Intel icc http://www.intel.com/
		CXX_VENDOR="Intel"
		CXXFLAGS="$CXXFLAGS -ansi" dnl -strict_ansi flag causes trouble
		CXX_OPTIMIZE_FLAGS="-O3 -Zp16 -ip -ansi_alias"
		CXX_DEBUG_FLAGS="-g -O0 -C"
		CXX_PROFIL_FLAGS="-pg"
	;;
	*g++*|*c++*)  dnl GNU C++  http://gcc.gnu.org/
		CXX_VENDOR="GNU" 
		GCC_V=`$CXX --version`
		gcc_version=`expr "$GCC_V" : '.* \(@<:@0-9@:>@\)\..*'`
		gcc_release=`expr "$GCC_V" : '.* @<:@0-9@:>@\.\(@<:@0-9@:>@\).*'`
		if test $gcc_version -lt "3" ; then
			CXXFLAGS="$CXXFLAGS -ftemplate-depth-40"
			CXX_OPTIMIZE_FLAGS="-O2 -funroll-loops -fstrict-aliasing -fno-gcse"
		else
			#CXXFLAGS=""
			CXX_OPTIMIZE_FLAGS="-O3 -funroll-loops -fstrict-aliasing -fomit-frame-pointer -ffast-math"
		fi
		CXX_DEBUG_FLAGS="-g"
		CXX_PROFIL_FLAGS="-pg"
	;;
	*pgCC*) dnl Portland Group   http://www.pgroup.com
		CXX_VENDOR="PGI"
		#CXXFLAGS=""
		CXX_OPTIMIZE_FLAGS="-O4 -Mnoframe -Mnodepchk -Minline=levels:25"
		CXX_DEBUG_FLAGS="-g -O0"
	;;
	*pathCC*) dnl Pathscale pathCC compiler   http://www.pathscale.com
		CXX_VENDOR="pathCC"
		CXXFLAGS="$CXXFLAGS -D__PATHSCALE -ansi"
		CXX_OPTIMIZE_FLAGS="-O3 -fstrict-aliasing -finline-functions"
		CXX_DEBUG_FLAGS="-g"
		CXX_PROFIL_FLAGS="-pg"
		#AR="$CXX"
		#AR_FLAGS="-ar -o"
	;;
	*) 
		ac_cxx_flags_preset=no
	;;
	esac
	AC_MSG_RESULT([yes])
else
	AC_MSG_RESULT([no])
fi

if test "$ac_cxx_flags_preset" = yes ; then
		AC_MSG_NOTICE([Setting compiler flags for $CXX_VENDOR $CXX])
else

	AC_MSG_NOTICE([No flags preset found for $CXX])
fi

])


