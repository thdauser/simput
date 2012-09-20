if FORTRAN
  FSRC=fsimput.f95
  FFLS=-Df2cFortran 
else
  FSRC=
  FFLS=
endif

AM_CFLAGS = -I@top_srcdir@/extlib/cfitsio -I@top_srcdir@/extlib/wcslib/C $(FFLS) 
#-fPIC

########### LIBRARIES #################

# Libraries that are installed in the 'lib' directory.
lib_LTLIBRARIES= libsimput.la

# Sources:
libsimput_la_SOURCES= datastruct.c fileaccess.c datahandling.c $(FSRC)

############ HEADERS #################

# Define the headers that will be installed in $(includedir):
include_HEADERS= simput.h headas.h headas_utils.h headas_stdio.h heasp.h
noinst_HEADERS = common.h
