dnl
dnl
dnl
AC_DEFUN([ARAN_C_VARIABLE_ARRAY], [
  AC_MSG_CHECKING([for C compiler variable length array support])
  AC_LANG_PUSH(C)dnl
  AC_COMPILE_IFELSE([
    AC_LANG_SOURCE([[ void _vlength (int n) {char buf[n];} ]])
  ] , [AC_MSG_RESULT(ok)
  ] , [AC_MSG_FAILURE([Variable length arrays are not supported by this compiler.])
  ])dnl
  AC_LANG_POP
])dnl
