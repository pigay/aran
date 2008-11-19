dnl
dnl Find the absolute path to C99 complex.h for use C++ programs that use
dnl LibAran
dnl
AC_DEFUN([ARAN_CXX_CCOMPLEX], [
  AC_MSG_NOTICE([Checking if C++ provides <ccomplex> header])
  cat > conftest.c << EOF
#include <complex.h>
EOF

  cat > conftest.sed << EOF
s/.*@<:@ @:>@\(@<:@/a-zA-Z0-9._@:>@\+complex\.h\).*/\1/p
d
EOF

  $CC $CFLAGS -M conftest.c > conftest.dep
  aran_cc_complex_h_path=`sed -f conftest.sed conftest.dep`
  rm -f conftest.*

  if test -r $aran_cc_complex_h_path ; then
    ARAN_CC_COMPLEX_H_INCLUDE="#include <$aran_cc_complex_h_path>"
  else
    AC_MSG_WARN([configure was not able to find absolute path to complex.h])
    ARAN_CC_COMPLEX_H_INCLUDE="#include <ccomplex>"
  fi

  AC_SUBST(ARAN_CC_COMPLEX_H_INCLUDE)
])
