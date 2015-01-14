# First, include the real Makefile
include Makefile

# Then, define the other targets needed by Automake Makefiles.

version=3.3.0 #`grep -R "*version = (float)*" *.c | awk '{print substr($$5,0,5);}'`
filename=libcfitsio.so
CFITSIO_PREFIX = cfitsio
CFITSIO_INCLUDE = ${DESTDIR}${prefix}/include/${CFITSIO_PREFIX}
INSTALL_DIRS_ALL    = ${INSTALL_DIRS} ${CFITSIO_INCLUDE}

inc_files=fitsio.h fitsio2.h longnam.h drvrsmem.h

install: ${filename} $(INSTALL_DIRS_ALL)
	make install-inc
	make install-lib

install-lib:
	/bin/cp ${filename} ${CFITSIO_LIB}/${filename}.${version}
	cd ${CFITSIO_LIB} && ln -fs ${filename}.${version} ${filename}

install-inc: $(INSTALL_DIRS_ALL)
	@if [ ! -d ${CFITSIO_INCLUDE} ]; then mkdir -p ${CFITSIO_INCLUDE}; fi
	/bin/cp ${inc_files}  ${CFITSIO_INCLUDE}/
	$(foreach inc,$(inc_files), cd ${CFITSIO_INCLUDE}/../ && ln -fs ${CFITSIO_PREFIX}/$(inc); )

$(INSTALL_DIRS_ALL):
	@if [ ! -d $@ ]; then mkdir -p $@; fi
