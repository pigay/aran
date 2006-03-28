dnl
dnl
dnl
AC_DEFUN([ARAN_C_COMPLEX], [
  AC_MSG_CHECKING([for C compiler complex numbers support])
  AC_LANG_PUSH(C)dnl
  AC_COMPILE_IFELSE([
    AC_LANG_PROGRAM([[ #include <complex.h> ]],
                    [[ float complex fc = I; double complex dc = I; ]])
  ] , [AC_MSG_RESULT(ok)
  ] , [AC_MSG_FAILURE([C complex numbers not supported.])
  ])dnl
  AC_LANG_POP
])dnl
