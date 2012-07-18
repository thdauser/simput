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

