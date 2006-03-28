# require automake-1.7
AUTOMAKE_OPTIONS = 1.7

SUBDIRS = gencode aran m4 tests docs

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = aran.pc

$(pkgconfig_DATA): config.status

EXTRA_DIST = autogen.sh types.list

ACLOCAL_AMFLAGS = -I m4
