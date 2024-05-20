# Copyright 2024 Tristan Britt
# distributed under the gpl_v3 license
# email: hello@tbritt.xyz
AC_INIT([cestimator], [0.1], [hello@tbritt.xyz], [cestimator])
AC_PREREQ(2.64)
AM_INIT_AUTOMAKE([foreign])

# check system type (no cross-compilation for now)
AC_CANONICAL_BUILD

# Add all needed -D options to dflags
dflags=""

AC_ARG_WITH(visualization,
   [AS_HELP_STRING([--with-visualization],
       [(no|yes|<path>) Use visualization, if yes configure assumes that a valid installation of MathGL is available. (default: no)])],
   [if  test "$withval" = "no" ; then
      with_mgl=0
   else
      with_mgl=1
   fi],
   [with_mgl=0])


AC_LANG_PUSH(C++)

if test "$with_mgl" -eq 1 ; then
    dflags="$dflags -D__VISUALIZER"
fi

AC_LANG_POP(C++)

AC_SUBST(dflags)

AC_CONFIG_FILES([Makefile])

AC_OUTPUT