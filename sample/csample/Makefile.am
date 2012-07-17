AM_CFLAGS = -I@top_srcdir@/libsimput 
AM_CFLAGS+= -Wall -Werror

############ BINARIES #################

# The following line lists the programs that should be created.
noinst_PROGRAMS = csample

csample_SOURCES = csample.c
csample_LDADD = @top_builddir@/libsimput/libsimput.la

########### Sample data ##############

sampledir = $(pkgdatadir)
dist_sample_DATA = image1.fits image2.fits
