AM_FCFLAGS = -I@top_srcdir@/libsimput 
AM_FCFLAGS+= -Wall

############ BINARIES #################

if FORTRAN
  PROGS=fsample
else
  PROGS=
endif

# The following line lists the programs that should be created.
noinst_PROGRAMS=$(PROGS)

fsample_SOURCES= fsample.f95
fsample_LDADD  = @top_builddir@/libsimput/libsimput.la 
fsample_LDADD += @top_builddir@/extlib/cfitsio/libcfitsio.a
fsample_LDADD += @top_builddir@/extlib/wcslib/C/libwcs-4.13.4.a
fsample_LDADD += @top_builddir@/extlib/fftw/libfftw3.la
