AM_CFLAGS = -I@top_srcdir@/libsimput 
AM_CFLAGS += -Wall 
#AM_CFLAGS += -Wall -Werror

############ BINARIES #################

# The following line lists the programs that should be created and stored
# in the 'bin' directory.
noinst_PROGRAMS = sample

sample_SOURCES = sample.c
sample_LDADD  = @top_builddir@/libsimput/libsimput.la



########### Sample data ##############

#sampledir = $(pkgdatadir)
#dist_sample_DATA = image1.fits image2.fits
