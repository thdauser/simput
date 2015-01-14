AM_CFLAGS =-I@top_srcdir@/extlib/cfitsio 
AM_CFLAGS+=-I@top_srcdir@/extlib/wcslib/C
AM_CFLAGS+=-I@top_srcdir@/extlib/ape/include
AM_CFLAGS+=-I@top_srcdir@/extlib/atFunctions/include
AM_CFLAGS+=-I@top_srcdir@/extlib/fftw/api

########## DIRECTORIES ###############

# Directory where to install the PIL parameter files.
pfilesdir=$(pkgdatadir)/pfiles
dist_pfiles_DATA=labnh.par

# Include README in the distributed package.
EXTRA_DIST=README

########### LIBRARIES ################

# Libraries that are installed in the 'lib' directory.
lib_LTLIBRARIES=liblabnh.la libposstring.la

# Sources:
liblabnh_la_SOURCES=liblabnh.c
libposstring_la_SOURCES=libposstring.c

liblabnh_la_LDFLAGS = -version-info 1:0:0
libposstring_la_LDFLAGS = -version-info 1:0:0

############ HEADERS #################

# Define the headers that will be installed in $(includedir):
include_HEADERS=labnh.h posstring.h

############ BINARIES #################

# The following line lists the programs that should be created and
# stored in the 'bin' directory.
bin_PROGRAMS=labnh

labnh_SOURCES=labnh.c labnh.h
labnh_LDADD  =liblabnh.la libposstring.la
labnh_LDADD +=@top_builddir@/extlib/ape/src/libape.la
labnh_LDADD +=@top_builddir@/extlib/atFunctions/src/libatFunctions.la
