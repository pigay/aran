
doc:

AM_CFLAGS = @CFLAGS@ -DG_LOG_DOMAIN=\"Aran\"

lib_LTLIBRARIES = libaran.la
bin_PROGRAMS = aranprofile

libaran_la_built_headers = arancoefficientbufferl.h \
arancoefficientbufferf.h arancoefficientbufferd.h \
arancoefficientbufferc64.h arancoefficientbufferc128.h \
aranbinomialbufferl.h \
aranbinomialbufferf.h aranbinomialbufferd.h \
aranbinomialbufferc64.h aranbinomialbufferc128.h \
aranccomplex.h

libaran_la_built_sources = arancoefficientbufferl.c \
arancoefficientbufferf.c arancoefficientbufferd.c \
arancoefficientbufferc64.c arancoefficientbufferc128.c \
aranbinomialbufferl.c \
aranbinomialbufferf.c aranbinomialbufferd.c \
aranbinomialbufferc64.c aranbinomialbufferc128.c arancomplex.c


libaran_la_sources = aran.c aransolver2d.c aranbinomial.c \
aranlaurentseriesd.c arandevelopment2d.c aranlegendre.c \
aransphericalharmonic.c aransphericalseriesd.c arandevelopment3d.c \
aransolver3d.c aranwigner.c aranwignerrepo.c aransphericalseriesd-translate.c \
aransphericalseriesd-kkylin.c aransphericalseriesd-rotate.c aranpoly1d.c \
aranlinear.c aranfit.c aranpolynomialfit.c aranprofile.c aranrusage.c \
aranprofiledb.c

libaran_la_headers = arancomplex.h aran.h aransolver2d.h aranbinomial.h \
aranlaurentseriesd.h arandevelopment2d.h aranlegendre.h \
aransphericalharmonic.h aransphericalseriesd.h arandevelopment3d.h \
aransolver3d.h aranwigner.h aranwignerrepo.h aranpoly1d.h aranlinear.h \
aranfit.h aranpolynomialfit.h aranprofile.h aranrusage.h aranprofiledb.h

libaran_la_noinst_headers = aransphericalseriesd-private.h aranwigner-private.h

libaran_la_SOURCES = $(libaran_la_built_headers) $(libaran_la_built_sources) \
$(libaran_la_headers) $(libaran_la_sources) $(libaran_la_noinst_headers)

aranincludedir = $(includedir)/aran
araninclude_HEADERS = $(libaran_la_built_headers) $(libaran_la_headers)

libaran_la_LIBADD = 
libaran_la_CPPFLAGS=-DARAN_DATA_DIR="\"$(pkgdatadir)\""

aranprofile_LDADD = libaran.la
aranprofile_SOURCES = aranprofile-main.c 

noinst_HEADERS = $(libaran_la_noinst_headers)

noinst_DATA = arancoefficientbuffer@t@.h arancoefficientbuffer@t@.c \
aranbinomialbuffer@t@.h aranbinomialbuffer@t@.c aranccomplex.h.in

pkgdata_DATA = profiledb.ini

EXTRA_DIST = $(noinst_DATA) $(pkgdata_DATA)
BUILT_SOURCES = $(libaran_la_built_headers) $(libaran_la_built_sources)

# rules for building profile db
profiledb.ini :
	./aranprofile > $@

all : aranprofile profiledb.ini

profile : clean_profile aranprofile profiledb.ini

clean_profile :
	\rm profiledb.ini

.PHONY: profile clean_profile

# rules for building generated sources
types_list = $(top_srcdir)/types.list

%l.h %d.h %f.h %c64.h %c128.h : %@t@.h ../gencode/typegens.py $(types_list)
	$(PYTHON) ../gencode/typegens.py -f=$(types_list) $<

%l.c %d.c %f.c %c64.c %c128.c : %@t@.c ../gencode/typegens.py $(types_list)
	$(PYTHON) ../gencode/typegens.py -f=$(types_list) $<

%d-private.h %f-private.h : %@t@-private.h ../gencode/typegens.py $(types_list)
	$(PYTHON) ../gencode/typegens.py -f=$(types_list) $<

