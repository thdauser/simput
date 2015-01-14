AM_CFLAGS =-I@top_srcdir@/extlib/cfitsio 

########### LIBRARIES #################

# Libraries that are installed in the 'lib' directory.
lib_LTLIBRARIES=libhdio.la

# Sources:
libhdio_la_SOURCES=headas_stdio.c
libhdio_la_LDFLAGS= -version-info 1:0:0

############ HEADERS #################

# Define the headers that will be installed in $(includedir):
include_HEADERS=headas_stdio.h
