
AM_CFLAGS = @CFLAGS@ -DG_LOG_DOMAIN=\"Aran-Tests\"

INCLUDES = -I$(top_srcdir) -I$(top_srcdir)/aran

noinst_PROGRAMS = 

noinst_PROGRAMS += dummypot taylor laurent legendre sphericalharmonic \
sphericalseriesd multipole3 taylor3 m2l3 newtonpot3 dev3 special_legendre \
sphericalharmonic_pregradient gradient3 newtonfield3 wigner rotation \
dummypotparallel newtonfield3parallel profiledb

LDADD = $(top_srcdir)/aran/libaran.la

$(srcdir)/package.m4: $(top_srcdir)/configure.ac
	{                                      \
	echo '# Signature of the current package.'; \
	echo 'm4_define([AT_PACKAGE_NAME],      [@PACKAGE_NAME@])'; \
	echo 'm4_define([AT_PACKAGE_TARNAME],   [@PACKAGE_TARNAME@])'; \
	echo 'm4_define([AT_PACKAGE_VERSION],   [@PACKAGE_VERSION@])'; \
	echo 'm4_define([AT_PACKAGE_STRING],    [@PACKAGE_STRING@])'; \
	echo 'm4_define([AT_PACKAGE_BUGREPORT], [@PACKAGE_BUGREPORT@])'; \
	} >$(srcdir)/package.m4

testscripts = dummypot.at taylor.at laurent.at legendre.at \
sphericalharmonic.at sphericalseriesd.at multipole3.at taylor3.at m2l3.at \
newtonpot3.at dev3.at special_legendre.at sphericalharmonic_pregradient.at \
gradient3.at newtonfield3.at wigner.at rotation.at profiledb.at

profiledbs = profiledb-dummypot.ini profiledb-newtonpot3.ini \
profiledb-newtonfield3.ini

EXTRA_DIST = testsuite.at $(testscripts) testsuite profiledb-test.ini \
$(profiledbs)

TESTSUITE = $(srcdir)/testsuite
check-local: atconfig $(TESTSUITE)
	$(SHELL) $(TESTSUITE)

AUTOTEST = $(AUTOM4TE) --language=autotest
$(TESTSUITE): $(srcdir)/testsuite.at $(srcdir)/package.m4 $(testscripts) 
	$(AUTOTEST) -I $(srcdir) $@.at -o $@.tmp
	mv $@.tmp $@

# rules for building profile db
profiledb-dummypot.ini : 
	./dummypot -profile > $@

profiledb-newtonpot3.ini : 
	./newtonpot3 -profile > $@

profiledb-newtonfield3.ini : 
	./newtonfield3 -profile > $@

all : dummypot $(profiledbs)

profile : clean_profile dummypot 
clean_profile :
	\rm $(profiledbs)

.PHONY: profile clean_profile

