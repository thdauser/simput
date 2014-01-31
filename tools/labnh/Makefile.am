AM_CFLAGS =-I@top_srcdir@/extlib/cfitsio 
AM_CFLAGS+=-I@top_srcdir@/extlib/wcslib/C
AM_CFLAGS+=-I@top_srcdir@/extlib/ape/include
AM_CFLAGS+=-I@top_srcdir@/extlib/atFunctions/include
AM_CFLAGS+=-I@top_srcdir@/extlib/fftw/api

########## DIRECTORIES ###############

# Directory where to install the PIL parameter files.
pfilesdir=@PFILESDIR@
dist_pfiles_DATA=labnh.par

# Include README in the distributed package.
EXTRA_DIST=README

########### LIBRARIES ################

# Libraries that are installed in the 'lib' directory.
lib_LIBRARIES=liblabnh.a libposstring.a

# Sources:
liblabnh_a_SOURCES=liblabnh.c
libposstring_a_SOURCES=libposstring.c

############ HEADERS #################

# Define the headers that will be installed in $(includedir):
include_HEADERS=labnh.h posstring.h

############ BINARIES #################

# The following line lists the programs that should be created and
# stored in the 'bin' directory.
bin_PROGRAMS=labnh

labnh_SOURCES=labnh.c labnh.h
labnh_LDADD  =liblabnh.a libposstring.a
labnh_LDADD +=@top_builddir@/extlib/ape/src/libape.a
labnh_LDADD +=@top_builddir@/extlib/atFunctions/src/libatFunctions.a
labnh_LDADD +=@top_builddir@/extlib/cfitsio/libcfitsio.a
labnh_LDADD +=@top_builddir@/extlib/wcslib/C/libwcs-4.13.4.a
