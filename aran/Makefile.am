
doc:

AM_CFLAGS = @CFLAGS@ -DG_LOG_DOMAIN=\"Aran\"

lib_LTLIBRARIES = libaran.la

libaran_la_built_headers = arancoefficientbufferl.h \
arancoefficientbufferf.h arancoefficientbufferd.h \
arancoefficientbufferc64.h arancoefficientbufferc128.h \
aranbinomialbufferl.h \
aranbinomialbufferf.h aranbinomialbufferd.h \
aranbinomialbufferc64.h aranbinomialbufferc128.h

libaran_la_built_sources = arancoefficientbufferl.c \
arancoefficientbufferf.c arancoefficientbufferd.c \
arancoefficientbufferc64.c arancoefficientbufferc128.c \
aranbinomialbufferl.c \
aranbinomialbufferf.c aranbinomialbufferd.c \
aranbinomialbufferc64.c aranbinomialbufferc128.c


libaran_la_sources = aran.c aransolver2d.c aranbinomial.c \
aranlaurentseriesd.c arandevelopment2d.c aranlegendre.c \
aransphericalharmonic.c aransphericalseriesd.c arandevelopment3d.c \
aransolver3d.c aranwigner.c aranwignerrepo.c aransphericalseriesd-translate.c \
aransphericalseriesd-kkylin.c aransphericalseriesd-rotate.c

libaran_la_headers = arancomplex.h aran.h aransolver2d.h aranbinomial.h \
aranlaurentseriesd.h arandevelopment2d.h aranlegendre.h \
aransphericalharmonic.h aransphericalseriesd.h arandevelopment3d.h \
aransolver3d.h aranwigner.h aranwignerrepo.h

libaran_la_noinst_headers = aransphericalseriesd-private.h aranwigner-private.h

libaran_la_SOURCES = $(libaran_la_built_headers) $(libaran_la_built_sources) \
$(libaran_la_headers) $(libaran_la_sources) $(libaran_la_noinst_headers)

aranincludedir = $(includedir)/aran
araninclude_HEADERS = $(libaran_la_built_headers) $(libaran_la_headers)

libaran_la_LIBADD = 

noinst_HEADERS = $(libaran_la_noinst_headers)

noinst_DATA = arancoefficientbuffer@t@.h arancoefficientbuffer@t@.c \
aranbinomialbuffer@t@.h aranbinomialbuffer@t@.c 

EXTRA_DIST = $(noinst_DATA)

BUILT_SOURCES = $(libaran_la_built_headers) $(libaran_la_built_sources)

types_list = $(top_srcdir)/types.list

%l.h %d.h %f.h %c64.h %c128.h : %@t@.h ../gencode/typegens.py $(types_list)
	$(PYTHON) ../gencode/typegens.py -f=$(types_list) $<

%l.c %d.c %f.c %c64.c %c128.c : %@t@.c ../gencode/typegens.py $(types_list)
	$(PYTHON) ../gencode/typegens.py -f=$(types_list) $<

%d-private.h %f-private.h : %@t@-private.h ../gencode/typegens.py $(types_list)
	$(PYTHON) ../gencode/typegens.py -f=$(types_list) $<

