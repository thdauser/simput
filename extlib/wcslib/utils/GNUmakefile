#-----------------------------------------------------------------------------
# GNU makefile for WCSLIB 4.13 utilities: fitshdr, HPXcvt, wcsgrid and wcsware.
#
# Summary of the main targets
# ---------------------------
#   build:     Build all utilities.
#   clean (or cleaner): Delete intermediate object files.
#   distclean (or realclean): cleaner, and also delete the executables.
#   cleanest: distclean, and also delete the man pages.
#
# Notes:
#   1: If you need to make changes then preferably modify ../makedefs.in
#      instead and re-run configure.
#
#   2: In compiling these utilities, this makefile assumes that the WCSLIB 4.13
#      sources reside in ../{pgsbox,C} (as in the distribution kit).
#
# Author: Mark Calabretta, Australia Telescope National Facility
# http://www.atnf.csiro.au/~mcalabre/index.html
# $Id: GNUmakefile,v 4.13.1.2 2012/03/21 02:58:52 cal103 Exp cal103 $
#-----------------------------------------------------------------------------
# Get configure settings.
include ../makedefs

UTILS := fitshdr

ifneq "$(CFITSIOINC)" ""
ifneq "$(CFITSIOLIB)" ""
  UTILS += HPXcvt wcsware

  ifneq "$(PGPLOTINC)" ""
  ifneq "$(PGPLOTLIB)" ""
    UTILS += wcsgrid
  endif
  endif
endif
endif

MAN := $(addsuffix .1,$(UTILS))

WCSLIB    := ../C/$(WCSLIB)
PGSBOXLIB := ../pgsbox/libpgsbox-$(LIBVER).a

ifneq "$(GETWCSTAB)" ""
  GETWCSTAB := ../C/$(GETWCSTAB)
endif

CPPFLAGS += -I.. -I../C

vpath %.h  ..:../C:../pgsbox
vpath %.in ..


# Static and static pattern rules
#--------------------------------

.PHONY : build clean cleaner cleanest distclean install man realclean

build : $(UTILS)

fitshdr : fitshdr.c
	-@ echo ''
	   $(CC) $(CPPFLAGS) $(CFLAGS) -c -o $@.o $<
	   $(LD) -o $@ $@.o $(LDFLAGS) $(LIBS)
	-@ $(RM) $@.o

HPXcvt : HPXcvt.c
	-@ echo ''
	   $(CC) $(CPPFLAGS) $(CFITSIOINC) $(CFLAGS) -c -o $@.o $<
	   $(LD) -o $@ $@.o $(LDFLAGS) $(CFITSIOLIB) $(LIBS)
	-@ $(RM) $@.o

wcsware : wcsware.c $(GETWCSTAB) $(WCSLIB)
	-@ echo ''
	   $(CC) -DDO_CFITSIO $(CPPFLAGS) -I../pgsbox -I../C $(CFITSIOINC) \
	     $(CFLAGS) -c -o $@.o $<
	   $(LD) -o $@ $@.o $(GETWCSTAB) $(LDFLAGS) $(CFITSIOLIB) $(WCSLIB) \
	     $(FLIBS) $(LIBS)
	-@ $(RM) $@.o

wcsgrid : wcsgrid.c $(GETWCSTAB) $(PGSBOXLIB) $(WCSLIB)
	-@ echo ''
	   $(CC) -DDO_CFITSIO $(CPPFLAGS) -I../pgsbox -I../C $(PGPLOTINC) \
	     $(CFITSIOINC) $(CFLAGS) -c -o $@.o $<
	   $(LD) -o $@ $@.o $(GETWCSTAB) $(LDFLAGS) $(PGSBOXLIB) \
	     $(PGPLOTLIB) $(CFITSIOLIB) $(WCSLIB) $(FLIBS) $(LIBS)
	-@ $(RM) $@.o

../C/getwcstab.o :
	   $(MAKE) -C ../C $(@F)

clean cleaner :
	-  $(RM) -r *.o *.i a.out core *.dSYM $(EXTRA_CLEAN)

distclean realclean : cleaner
	-  $(RM) $(UTILS)

cleanest : distclean
	-  $(RM) $(MAN)

$(PGSBOXLIB) :
	-@ echo ''
	   $(MAKE) -C ../pgsbox lib

$(WCSLIB) :
	-@ echo ''
	   $(MAKE) -C ../C lib

install : build
	-  if [ ! -d "$(BINDIR)" ] ; then \
	     $(INSTALL) -d -m 2775 $(BINDIR) ; \
	   fi
	   $(INSTALL) -m 755 $(UTILS) $(BINDIR)
	-  if [ ! -d "$(MANDIR)" ] ; then \
	     $(INSTALL) -d -m 2775 $(MANDIR)/man1 ; \
	   fi
	   $(INSTALL) -m 644 $(MAN) $(MANDIR)/man1

GNUmakefile : ../makedefs ;

../makedefs ../wcsconfig.h ../wcsconfig_utils.h : makedefs.in wcsconfig.h.in \
    wcsconfig_utils.h.in ../config.status
	-@ $(RM) ../wcsconfig.h ../wcsconfig_utils.h
	   cd .. && ./config.status

man : $(MAN)

fitshdr.1 : fitshdr
	help2man --no-discard-stderr --version-string=$(LIBVER) \
	  -n "List headers from a FITS file" -N ./$< > $@

wcsware.1 : wcsware
	help2man --no-discard-stderr --version-string=$(LIBVER) \
	  -n "Extract WCS keywords for an image" -N ./$< > $@

wcsgrid.1 : wcsgrid
	help2man --no-discard-stderr --version-string=$(LIBVER) \
	  -n "Extract WCS keywords for an image" -N ./$< > $@

HPXcvt.1 : HPXcvt
	help2man --no-discard-stderr --version-string=$(LIBVER) \
	  -n "Reorganise HEALPix data into a 2-D FITS image" -N ./$< > $@

# Dependency lists.
fitshdr : wcsconfig.h wcsconfig_utils.h
wcsware : getwcstab.h wcs.h wcsfix.h wcshdr.h
wcsgrid : cpgsbox.h getwcstab.h wcs.h wcsfix.h wcshdr.h
