# The sub-directories are built before the current directory.
# In order to change this, include "." in the list of SUBDIRS.
SUBDIRS=. ape atFunctions heainit heaio heautils heasp

DIST_SUBDIRS=ape atFunctions heainit heaio heautils heasp fftw
EXTRA_DIST=cfitsio wcslib 

# We assume that the configure-script is executed by the main configure.ac in SIMPUT

all-local:
	make all-wcslib
	make all-cfitsio
	make all-fftw

install-data-local: 
	make install-wcslib
	make install-cfitsio
	make install-fftw
	make install-headas-main

clean-local:
	make clean-fftw
	make clean-cfitsio
	make clean-wcslib


##### HEADAS MAIN #####
headas_main_dir=heainit
install-headas-main:
	/bin/cp ${headas_main_dir}/headas_main.c ${includedir}

##### CFITSIO #####

cfitsio_dir=cfitsio
all-cfitsio:
	cd ${cfitsio_dir} &&  make shared

# need to create and linke an adjusted Makefile
install-cfitsio:
	cd ${cfitsio_dir} && make install

clean-cfitsio:
	cd ${cfitsio_dir} && make clean

##### FFTW #####

fftw_dir=fftw
all-fftw:
	cd ${fftw_dir} && make 

install-fftw:
	cd ${fftw_dir} && make install

clean-fftw:
	cd ${fftw_dir} && make clean

##### WCS #####

wcslib_dir=wcslib
LIBVER=4.25.1

LIBST_OSX=dylib
LIBST_LINUX=so

if IS_OSX
SHRLIB=libwcs.${LIBVER}.${LIBST_OSX}
SHRLN=libwcs.${LIBST_OSX}
else
SHRLIB=libwcs.${LIBST_LINUX}.${LIBVER}
SHRLN=libwcs.${LIBST_LINUX}
endif

all-wcslib:
	cd ${wcslib_dir} && make && cd ./C/ && ln -sf ${SHRLIB} ${SHRLN}

INSTALL_DIRS = ${libdir}
install-wcslib: $(INSTALL_DIRS)
	cd ${wcslib_dir} && make install
	make ${wcslib_inc_files}

clean-wcslib:
	cd ${wcslib_dir} && make distclean

# just to create the links
wcslib_inc_files=cel.h      getwcstab.h  log.h  spc.h  spx.h  wcsconfig_f77.h  wcserr.h  wcs.h  \
		 wcslib.h   wcsprintf.h  wcsunits.h fitshdr.h  lin.h        prj.h  sph.h  tab.h \
		 wcsconfig.h      wcsfix.h  wcshdr.h  wcsmath.h  wcstrig.h    wcsutil.h
.PHONY: ${wcslib_inc_files}
${wcslib_inc_files}: 
	cd ${includedir} && ln -fs wcslib/$@			

$(INSTALL_DIRS):
	@if [ ! -d $@ ]; then mkdir -p $@; fi
