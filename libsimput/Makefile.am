#AM_CFLAGS  = -I${HEADAS}/include
#AM_CFLAGS += -Wall 
#AM_CFLAGS += -Wall -Werror 


########## DIRECTORIES ###############

# Directory where to install the FITS template files.
#fitstemplatesdir = $(pkgdatadir)/templates

# Directory where to install the XML files.
#xmldir = $(pkgdatadir)/xml

########### LIBRARIES #################

# Libraries that are installed in the 'lib' directory.
lib_LTLIBRARIES = libsimput.la

# Sources:
libsimput_la_SOURCES = simput.c

############ HEADERS #################

# Define the headers that will be installed in $(includedir):
include_HEADERS = simput.h






############ BINARIES #################

# The following line lists the programs that should be created.
noinst_PROGRAMS = simput_fortran

simput_fortran_SOURCES = simput_fortran.f95
simput_fortran_LDADD  = @top_builddir@/libsimput/libsimput.la

