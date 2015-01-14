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
dist_pfiles_DATA=simputpsd.par

############ BINARIES #################

# The following line lists the programs that should be created and
# stored in the 'bin' directory.
bin_PROGRAMS=simputpsd

simputpsd_SOURCES=simputpsd.c simputpsd.h
simputpsd_LDADD =@top_builddir@/libsimput/libsimput.la
simputpsd_LDADD+=@top_builddir@/extlib/heainit/libhdinit.la
simputpsd_LDADD+=@top_builddir@/extlib/heaio/libhdio.la
simputpsd_LDADD+=@top_builddir@/extlib/heautils/libhdutils.la
simputpsd_LDADD+=@top_builddir@/extlib/heasp/libhdsp.la
simputpsd_LDADD+=@top_builddir@/extlib/ape/src/libape.la
