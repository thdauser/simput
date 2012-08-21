AM_CFLAGS = -I@top_srcdir@/libsimput -I@top_srcdir@/extlib/cfitsio -I@top_srcdir@/extlib/wcslib/C
AM_CFLAGS+= -Wall -Werror

############ BINARIES #################

# The following line lists the programs that should be created.
noinst_PROGRAMS= csample

csample_SOURCES= csample.c
csample_LDADD  = @top_builddir@/libsimput/libsimput.la 
#csample_LDADD += @top_builddir@/extlib/cfitsio/libcfitsio.a 
csample_LDADD += @top_builddir@/extlib/wcslib/C/libwcs-4.13.4.a
csample_LDADD += @top_builddir@/extlib/fftw/libfftw3.la

########### Sample data ##############

sampledir       = $(pkgdatadir)
dist_sample_DATA= image1.fits image2.fits
