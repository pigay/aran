# Process this file with autoconf to produce a configure script.

# glib style versioning
#
# Making releases:
#   aran_micro_version += 1;
#   aran_interface_age += 1;
#   aran_binary_age += 1;
# if any functions have been added, set aran_interface_age to 0.
# if backwards compatibility has been broken,
# set aran_binary_age _and_ aran_interface_age to 0.
#
m4_define([aran_major_version], [0])
m4_define([aran_minor_version], [3])
m4_define([aran_micro_version], [0])
m4_define([aran_interface_age], [0])
m4_define([aran_binary_age],
          [m4_eval(100 * aran_minor_version + aran_micro_version)])
m4_define([aran_version],
          [aran_major_version.aran_minor_version.aran_micro_version])

# libtool version related macros
m4_define([aran_lt_release], [aran_major_version.aran_minor_version])
m4_define([aran_lt_current],
          [m4_eval(100 * aran_minor_version + aran_micro_version - aran_interface_age)])
m4_define([aran_lt_revision], [aran_interface_age])
m4_define([aran_lt_age], [m4_eval(aran_binary_age - aran_interface_age)])
m4_define([aran_lt_current_minus_age],
          [m4_eval(aran_lt_current - aran_lt_age)])


AC_INIT(aran, [aran_version])

AC_CONFIG_SRCDIR([aran/aran.c])

ARAN_MAJOR_VERSION=aran_major_version
ARAN_MINOR_VERSION=aran_minor_version
ARAN_MICRO_VERSION=aran_micro_version
ARAN_INTERFACE_AGE=aran_interface_age
ARAN_BINARY_AGE=aran_binary_age
ARAN_VERSION=aran_version

AC_SUBST(ARAN_MAJOR_VERSION)
AC_SUBST(ARAN_MINOR_VERSION)
AC_SUBST(ARAN_MICRO_VERSION)
AC_SUBST(ARAN_VERSION)
AC_SUBST(ARAN_INTERFACE_AGE)
AC_SUBST(ARAN_BINARY_AGE)

AC_DEFINE(ARAN_MAJOR_VERSION, [aran_major_version],
	  [Define to the ARAN major version])
AC_DEFINE(ARAN_MINOR_VERSION, [aran_minor_version],
	  [Define to the ARAN minor version])
AC_DEFINE(ARAN_MICRO_VERSION, [aran_micro_version],
	  [Define to the ARAN micro version])
AC_DEFINE(ARAN_INTERFACE_AGE, [aran_interface_age],
	  [Define to the ARAN interface age])
AC_DEFINE(ARAN_BINARY_AGE, [aran_binary_age],
	  [Define to the ARAN binary age])

# libtool versioning
LT_RELEASE=aran_lt_release
LT_CURRENT=aran_lt_current
LT_REVISION=aran_lt_revision
LT_AGE=aran_lt_age
LT_CURRENT_MINUS_AGE=aran_lt_current_minus_age
AC_SUBST(LT_RELEASE)
AC_SUBST(LT_CURRENT)
AC_SUBST(LT_REVISION)
AC_SUBST(LT_AGE)
AC_SUBST(LT_CURRENT_MINUS_AGE)

AC_CANONICAL_HOST
AC_CANONICAL_TARGET

AM_INIT_AUTOMAKE(aran, aran_version, [no-define])

AM_CONFIG_HEADER(aran-config.h)

# Checks for programs.
AC_PROG_CC
AC_PROG_AWK
AC_PROG_CC
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_MAKE_SET
AC_ISC_POSIX

AM_PROG_LEX
AC_PROG_YACC

AM_PROG_LIBTOOL

AM_MISSING_PROG([AUTOM4TE], [autom4te])

AC_CHECK_PROG(INDENT, indent, indent)

AM_PATH_PYTHON(2.2)

VL_PROG_CC_WARNINGS

# language choice
AC_LANG_C

AC_ARG_ENABLE(debug,
              [AC_HELP_STRING([--enable-debug], [turns on debugging])],
              [enable_debug=yes])

if test "x$enable_debug" = "xyes"; then
  case " $CFLAGS " in
    *[[\ \	]]-g[[\ \	]]*) ;;
    *) CFLAGS="$CFLAGS -g" ;;
  esac
fi

# Checks for libraries.
AC_CHECK_LIB(m, cos)


# Check for variable length array support
ARAN_C_VARIABLE_ARRAY

# We need libmingwex.a for complex numbers support on cygwin
case $host in
  *-*-mingw32* | *-*-cygwin*)
    mingw_include="-I/usr/include/mingw"
    mingw_ldflags="-L/usr/lib/mingw"
    mingw_libadd="-lcrtdll -lmingwex"
    CFLAGS="$CFLAGS $mingw_include"
    CXXFLAGS="$CXXFLAGS $mingw_include"
    LDFLAGS="$LDLFAGS $mingw_ldflags"
    AC_CHECK_LIB([crtdll], [_fopen])
    AC_CHECK_LIB([mingwex], [cexp])
    ;;
  *)
    AC_CHECK_LIB([m],[cexp])
    ;;
esac
# subst for pkg-config files
AC_SUBST(mingw_include)
AC_SUBST(mingw_ldflags)
AC_SUBST(mingw_libadd)

# Check for CC complex numbers support
ARAN_C_COMPLEX

# Check for C complex numbers support with C++
ARAN_CXX_CCOMPLEX

PKG_CHECK_MODULES(BASE_DEPENDENCIES, [glib-2.0 >= 2.0.0 gobject-2.0 >= 2.0.0
                                      gmodule-2.0 >= 2.0.0 vsgd >= 0.2.0])

CFLAGS="$BASE_DEPENDENCIES_CFLAGS $CFLAGS"
LIBS="$BASE_DEPENDENCIES_LIBS $LIBS "
AC_SUBST(CFLAGS)
AC_SUBST(LIBS)

##################################################
# Check for gtk-doc.
##################################################

AC_ARG_WITH(html-dir, [  --with-html-dir=PATH path to installed docs ])

if test "x$with_html_dir" = "x" ; then
  HTML_DIR='${datadir}/gtk-doc/html'
else
  HTML_DIR=$with_html_dir
fi

AC_SUBST(HTML_DIR)

AC_CHECK_PROG(GTKDOC, gtkdoc-mkdb, true, false)

gtk_doc_min_version=0.6
if $GTKDOC ; then 
    gtk_doc_version=`gtkdoc-mkdb --version`
    AC_MSG_CHECKING([gtk-doc version ($gtk_doc_version) >= $gtk_doc_min_version])
    if perl <<EOF ; then
      exit (("$gtk_doc_version" =~ /^[[0-9]]+\.[[0-9]]+$/) &&
            ("$gtk_doc_version" >= "$gtk_doc_min_version") ? 0 : 1);
EOF
      AC_MSG_RESULT(yes)
   else
      AC_MSG_RESULT(no)
      GTKDOC=false
   fi
fi

# Let people disable the gtk-doc stuff.
AC_ARG_ENABLE(gtk-doc, [  --enable-gtk-doc  Use gtk-doc to build documentation [default=auto]], enable_gtk_doc="$enableval", enable_gtk_doc=auto)

if test x$enable_gtk_doc = xauto ; then
  if test x$GTKDOC = xtrue ; then
    enable_gtk_doc=yes
  else
    enable_gtk_doc=no 
  fi
fi

AM_CONDITIONAL(ENABLE_GTK_DOC, test x$enable_gtk_doc = xyes)

# Checks for header files.
AC_HEADER_STDC
# 


# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_TYPE_SIZE_T
AC_C_INLINE

# Checks for library functions.


# Testsuite configuration.
AC_CONFIG_TESTDIR(tests)


# Generate files.
AC_OUTPUT([
aran.pc
Makefile
gencode/Makefile
aran/Makefile
aran/aranccomplex.h
m4/Makefile
tests/Makefile
docs/Makefile
docs/reference/Makefile
docs/reference/version.xml
] )
