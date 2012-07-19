########### LIBRARIES #################

# Libraries that are installed in the 'lib' directory.
lib_LTLIBRARIES= libsimput.la

if FORTRAN
  FSRC=fsimput.f95
else
  FSRC=
endif

# Sources:
libsimput_la_SOURCES= simput.c $(FSRC)

############ HEADERS #################

# Define the headers that will be installed in $(includedir):
include_HEADERS= simput.h


