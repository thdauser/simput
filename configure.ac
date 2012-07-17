#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.65)
AC_INIT([simput],[0.22],[christian.schmid@sternwarte.uni-erlangen.de])
AC_CONFIG_AUX_DIR([build-aux]) # Store auxiliary files in the sub-directory 'build-aux'.
AM_INIT_AUTOMAKE([foreign -Wall -Werror])
AC_PROG_CC # Search for a C compiler and store it in the Variable 'CC'.
AC_PROG_FC # Search for a Fortran compiler.
AC_CONFIG_SRCDIR([libsimput/simput.c])

LT_INIT #([disable-shared]) #([shared]) 

# Compiler and Linker Flags:
# "-g -W -Wall -O0" for the compiler for debugging
# "-g" to create debugging information (seems to be used automatically by autotools)
# "-Werror" for the compiler to treat warnings as errors
# "-pg" for the compiler and linker to include information for gprof (Profiler)
# "-std=c99 -D_GNU_SOURCE" for C99 syntax standards
# " -O0" for zero optimization
CFLAGS="$CFLAGS -g -I$HEADAS/include -std=c99 -D_GNU_SOURCE -W -Wall -Df2cFortran"
# -DHEASP_CPP (required for HEADAS versions >=6.12)
LDFLAGS="$LDFLAGS -g -L$HEADAS/lib -W -Wall"
# FCFLAGS="-O0"

# Check for libraries:  (objdump -R/-T libncurses...)
AC_SEARCH_LIBS(sin, m, [], [AC_MSG_ERROR([No math library found!])], [])

####################################
# cfitsio
AC_SEARCH_LIBS(ffopen, 
	[cfitsio_3.29 cfitsio_3.27 cfitsio_3.25 cfitsio_3.24 cfitsio_3.22 cfitsio_3.181 cfitsio_3.12 cfitsio_3.07 cfitsio_3.04 cfitsio], 
	[], [AC_MSG_ERROR([No supported version of cfitsio library found!])], [])
#lib_cfitsio="no"
#AC_MSG_CHECKING([for versioned libcfitsio in \$HEADAS])
#my_cfitsio=`ls "${HEADAS}"/lib/libcfitsio*.so 2> /dev/null`
#if ! test -z "$my_cfitsio"
#then
#   my_cfitsio=`basename "$my_cfitsio"`
#   my_cfitsio=`expr "$my_cfitsio" : 'lib\(cfitsio_.*\)\.so'`
#fi
#if ! test -z "$my_cfitsio"
#then
#   AC_MSG_RESULT([yes ($my_cfitsio)])
#else
#   AC_MSG_RESULT([no])
#fi
#AC_SEARCH_LIBS(ffopen, cfitsio $my_cfitsio, [],
#          [AC_MSG_ERROR([No cfitsio library found!])], [])
##AC_SEARCH_LIBS(ffopen, cfitsio $my_cfitsio,
##          [lib_cfitsio=-l$ac_lib],
##          [AC_MSG_ERROR([No cfitsio library found!])],
##          [-lm])
#AC_SUBST(lib_cfitsio)

#################################### 
# WCS library:
AC_SEARCH_LIBS(wcspih, [wcs_6.11 wcs_6.10 wcs_6.9 wcs_6.8 wcs_6.7 wcs_6.6 \
		        wcs-4.4.4 wcs-4.3.3 wcs-4.3], [], \
		       [AC_MSG_ERROR([No WCS library found!])], [])
#################################### 
# FFTW library:
AC_SEARCH_LIBS(fftw_execute, [fftw3], [], \
		      	 [AC_MSG_ERROR([No fftw3 library found!])], [])
#################################### 
# Check for header files:
AC_HEADER_STDC
AC_CHECK_HEADERS([limit.h math.h stdio.h stdlib.h stdarg.h string.h])
#AC_CHECK_HEADERS([headas.h headas_error.h]) 
# -> problem: not accepted by preprocessor (but by compiler)!

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST

# Checks for library functions.


# The following line tells 'configure' to create a header file 'config.h'
# that is filled with definitions from this script.
# The 'config.h' contains for example the name of the package and can be
# included in the programm sources.
AC_CONFIG_HEADER([config.h])

# Create the following Makefiles (including subdirectories).
# A Makefile that is not listed here will not be created by 'automake'!
AC_CONFIG_FILES([Makefile
		libsimput/Makefile
		sample/Makefile
		sample/csample/Makefile
		sample/fsample/Makefile])
AC_OUTPUT

