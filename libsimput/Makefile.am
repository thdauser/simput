if FORTRAN
  FSRC=fsimput.f95
  FFLS=-Df2cFortran 
else
  FSRC=
  FFLS=
endif

AM_CFLAGS =-I@top_srcdir@/extlib/heainit
AM_CFLAGS+=-I@top_srcdir@/extlib/heaio
AM_CFLAGS+=-I@top_srcdir@/extlib/heautils
AM_CFLAGS+=-I@top_srcdir@/extlib/heasp
AM_CFLAGS+=-I@top_srcdir@/extlib/cfitsio
AM_CFLAGS+=-I@top_srcdir@/extlib/wcslib/C $(FFLS)
AM_CFLAGS+=-I@top_srcdir@/extlib/fftw/api
AM_CFLAGS+=-I@top_srcdir@/extlib/ape/include

########### LIBRARIES #################

# Libraries that are installed in the 'lib' directory.
lib_LTLIBRARIES=libsimput.la

# Sources:
libsimput_la_SOURCES=datastruct.c fileaccess.c datahandling.c vector.c	\
                    arf.c rmf.c input.c $(FSRC)
libsimput_la_LIBADD=@top_builddir@/extlib/heasp/libhdsp.la

libsimput_la_LDFLAGS = -version-info 2:0:1

############ HEADERS #################

# Define the headers that will be installed in $(includedir):
include_HEADERS=simput.h common.h vector.h arf.h rmf.h input.h
