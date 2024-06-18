# Copyright 2024 Tristan Britt
# distributed under the gpl_v3 license
# email: hello@tbritt.xyz
AC_INIT([cestimator], [0.1], [hello@tbritt.xyz], [cestimator])
AC_PREREQ(2.64)
AC_CONFIG_AUX_DIR([build_aux])
AM_INIT_AUTOMAKE([foreign])

# check system type (no cross-compilation for now)
AC_CANONICAL_BUILD

# Add all needed -D options to dflags
dflags=""
# Add libraries if visualization
ldflags=""

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
AC_CHECK_HEADERS([eigen3/Eigen/Dense],
        [has_eigen=yes; break;])
AS_IF([test "x$has_eigen" != "xyes"],
   [AC_MSG_ERROR([Unable to find the Eigen headers])])

if test "$with_mgl" -eq 1 ; then
    AC_MSG_NOTICE([Building with visualization libraries])

    if test ! -d "external"; then 
      mkdir external
    fi 

    cd external 
    if test ! -d "matplotplusplus"; then 
      git clone https://github.com/alandefreitas/matplotplusplus.git
    fi

    cd matplotplusplus
    mkdir build 
    cd build
    cmake -DCMAKE_INSTALL_PREFIX:PATH="$(pwd)" -DMATPLOTPP_BUILD_EXAMPLES=OFF      \
    -DMATPLOTPP_BUILD_SHARED_LIBS=ON    \
    -DMATPLOTPP_BUILD_TESTS=OFF         \
    -DCMAKE_INTERPROCEDURAL_OPTIMIZATION=ON ..
    make -j 6
    make install
    cd ../../..
    dflags="$dflags -D__VISUALIZER"
    ldflags="$ldflags -fsanitize=undefined -L./external/matplotplusplus/build/lib -lmatplot -L./external/matplotplusplus/build/lib/Matplot++ -lnodesoup -lstdc++fs -lpng -ljpeg -ltiff -lpthread"
    iflags="$iflags -I./external/matplotplusplus/build/include"
fi

AC_LANG_POP(C++)

AC_SUBST(dflags)
AC_SUBST(ldflags)
AC_SUBST(iflags)

AC_CONFIG_FILES([Makefile])

AC_OUTPUT