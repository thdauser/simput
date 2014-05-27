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

########### LIBRARIES #################

# Libraries that are installed in the 'lib' directory.
lib_LIBRARIES=libsimput.a

# Sources:
libsimput_a_SOURCES=datastruct.c fileaccess.c datahandling.c vector.c	\
                    arf.c rmf.c $(FSRC)

############ HEADERS #################

# Define the headers that will be installed in $(includedir):
include_HEADERS=simput.h
noinst_HEADERS =common.h vector.h arf.h rmf.h
