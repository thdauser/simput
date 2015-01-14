AM_CFLAGS =-I@top_srcdir@/libsimput
AM_CFLAGS+=-I@top_srcdir@/labnh
AM_CFLAGS+=-I@top_srcdir@/extlib/cfitsio 
AM_CFLAGS+=-I@top_srcdir@/extlib/wcslib/C
AM_CFLAGS+=-I@top_srcdir@/extlib/ape/include
AM_CFLAGS+=-I@top_srcdir@/extlib/atFunctions/include
AM_CFLAGS+=-I@top_srcdir@/extlib/fftw/api
AM_CFLAGS+=-I${HEADAS}/include
AM_LDFLAGS =-L${HEADAS}/lib 

########## DIRECTORIES ###############

# Directory where to install the PIL parameter files.
pfilesdir=$(pkgdatadir)/pfiles

############ BINARIES #################

if GALABS
  PROGS=galabs
  dist_pfiles_DATA=galabs.par
else
  PROGS=
endif

bin_PROGRAMS=$(PROGS)

galabs_SOURCES=galabs.c spectree.c spectree.h
galabs_LDADD  =@top_builddir@/labnh/liblabnh.la 
galabs_LDADD +=@top_builddir@/labnh/libposstring.la
galabs_LDADD +=@top_builddir@/libsimput/libsimput.la 
galabs_LDADD +=-lXSFunctions 
galabs_LDADD +=-lXSModel 
galabs_LDADD +=-lXSUtil 
galabs_LDADD +=-lXS 
galabs_LDADD +=-lCCfits_2.4 
galabs_LDADD +=-lgfortran
galabs_LDADD +=@top_builddir@/extlib/ape/src/libape.la
galabs_LDADD +=@top_builddir@/extlib/atFunctions/src/libatFunctions.la
