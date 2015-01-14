AM_CFLAGS =-I@top_srcdir@/libsimput
AM_CFLAGS+=-I@top_srcdir@/extlib/heainit
AM_CFLAGS+=-I@top_srcdir@/extlib/heaio
AM_CFLAGS+=-I@top_srcdir@/extlib/heautils
AM_CFLAGS+=-I@top_srcdir@/extlib/heasp
AM_CFLAGS+=-I@top_srcdir@/extlib/cfitsio
AM_CFLAGS+=-I@top_srcdir@/extlib/wcslib/C
AM_CFLAGS+=-Wall -Werror

############ BINARIES #################

# The following line lists the programs that should be created.
noinst_PROGRAMS= csample

csample_SOURCES=csample.c
csample_LDADD =@top_builddir@/libsimput/libsimput.la
csample_LDADD+=@top_builddir@/extlib/heainit/libhdinit.la
csample_LDADD+=@top_builddir@/extlib/heaio/libhdio.la
csample_LDADD+=@top_builddir@/extlib/heautils/libhdutils.la
csample_LDADD+=@top_builddir@/extlib/heasp/libhdsp.la
csample_LDADD+=@top_builddir@/extlib/ape/src/libape.la

########### Sample data ##############

sampledir       =$(pkgdatadir)
dist_sample_DATA=image1.fits image2.fits
