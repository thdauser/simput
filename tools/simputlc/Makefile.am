AM_CFLAGS =-I@top_srcdir@/libsimput 
AM_CFLAGS+=-I@top_srcdir@/extlib/cfitsio 
AM_CFLAGS+=-I@top_srcdir@/extlib/wcslib/C
AM_CFLAGS+=-I@top_srcdir@/extlib/ape/include
AM_CFLAGS+=-I@top_srcdir@/extlib/heainit
AM_CFLAGS+=-I@top_srcdir@/extlib/heaio
AM_CFLAGS+=-I@top_srcdir@/extlib/heautils
AM_CFLAGS+=-I@top_srcdir@/extlib/heasp
AM_CFLAGS+=-I@top_srcdir@/extlib/fftw/api
AM_CFLAGS+=-Wall 

########## DIRECTORIES ###############

# Directory where to install the PIL parameter files.
pfilesdir=$(pkgdatadir)/pfiles
dist_pfiles_DATA=simputlc.par

############ BINARIES #################

# The following line lists the programs that should be created and
# stored in the 'bin' directory.
bin_PROGRAMS=simputlc

simputlc_SOURCES=simputlc.c simputlc.h
simputlc_LDADD =@top_builddir@/libsimput/libsimput.a
simputlc_LDADD+=@top_builddir@/extlib/heainit/libhdinit.a
simputlc_LDADD+=@top_builddir@/extlib/heaio/libhdio.a
simputlc_LDADD+=@top_builddir@/extlib/heautils/libhdutils.a
simputlc_LDADD+=@top_builddir@/extlib/heasp/libhdsp.a
simputlc_LDADD+=@top_builddir@/extlib/ape/src/libape.a
simputlc_LDADD+=@top_builddir@/extlib/wcslib/C/libwcs-4.13.4.a
simputlc_LDADD+=@top_builddir@/extlib/fftw/.libs/libfftw3.a
simputlc_LDADD+=@top_builddir@/extlib/cfitsio/libcfitsio.a
