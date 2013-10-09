#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.65)
AC_INIT([simput],[0.60],[christian.schmid@sternwarte.uni-erlangen.de])
# Store auxiliary files in the sub-directory 'build-aux'.
AC_CONFIG_AUX_DIR([build-aux]) 
AM_INIT_AUTOMAKE([foreign -Wall -Werror])
AC_PROG_CC # Search for a C compiler and store it in the Variable 'CC'.
AC_PROG_FC # Search for a Fortran compiler.
AC_PROG_RANLIB
AC_CONFIG_SRCDIR([libsimput/fileaccess.c])

# Compiler and Linker Flags:
# "-g -W -Wall -O0" for the compiler for debugging
# "-g" to create debugging information (seems to be used automatically by autotools)
# "-Werror" for the compiler to treat warnings as errors
# "-pg" for the compiler and linker to include information for gprof (Profiler)
# "-std=c99 -D_GNU_SOURCE" for C99 syntax standards
# "-O0" for zero optimization
CFLAGS="$CFLAGS -g -std=c99 -D_GNU_SOURCE -W -Wall -fPIC -Dg77Fortran"
LDFLAGS="$LDFLAGS -g -W -Wall"
# FCFLAGS="-O0"

# Conditional compilation of Fortran code.
AC_ARG_ENABLE([fortran], 
              [AS_HELP_STRING([--enable-fortran], 
	                      [enable compilation of Fortran code])],
	          [fortran=true], [fortran=false])
AM_CONDITIONAL([FORTRAN], [test x$fortran = xtrue])

# Conditional compilation of program 'galabs' (requires HEAdas).
AC_ARG_ENABLE([galabs], 
              [AS_HELP_STRING([--enable-galabs], 
	                      [enable compilation of binary 'galabs'])],
	          [galabs=true], [galabs=false])
AM_CONDITIONAL([GALABS], [test x$galabs = xtrue])

# Check for libraries.
AC_SEARCH_LIBS(sin, m, [], [AC_MSG_ERROR([library math not found!])], [])
# ncurses and readline (required by ape library).
AC_SEARCH_LIBS(beep, ncurses, [], [AC_MSG_ERROR([library ncurses not found!])], [])
AC_SEARCH_LIBS(readline, readline, [], [AC_MSG_ERROR([library readline not found!])], [])

#if test -z "$GALABS"
#then
# 	AC_SEARCH_LIBS(_gfortran_unlink, gfortran, [], [AC_MSG_ERROR([library gfortran not found!])], [])
# 	AC_SEARCH_LIBS(frame_dummy, CCfits_2.4, [], [AC_MSG_ERROR([library CCfits_2.4 not found!])], [])
# 	AC_SEARCH_LIBS(xgetin_, XS, [], [AC_MSG_ERROR([library XS not found!])], [])
# 	AC_SEARCH_LIBS(xs_read, XSUtil, [], [AC_MSG_ERROR([library XUtil not found!])], [])
# 	AC_SEARCH_LIBS(xsatbl, XSModel, [], [AC_MSG_ERROR([library XSModel not found!])], [])
#	AC_SEARCH_LIBS(xsphab_, XSFunctions, [], [AC_MSG_ERROR([library XSFunctions not found!])], [])
#fi

# Check for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([limit.h math.h stdio.h stdlib.h stdarg.h string.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST

# PIL parameter files have to be installed to the $PFILES directory.
AC_ARG_VAR([PFILESDIR], [target directory for .par files])
if test -z "$PFILESDIR"
then
	if test -z "$HEADAS"
	then
		AC_MSG_ERROR([Environment variable '\$HEADAS' is not set], [])
	else
		PFILESDIR="\${HEADAS}/pfiles"
	fi
fi
AC_SUBST([PFILESDIR])

# The following line tells 'configure' to create a header file 'config.h'
# that is filled with definitions from this script.
# The 'config.h' contains for example the name of the package and can be
# included in the programm sources.
AC_CONFIG_HEADER([config.h])

# Create the following Makefiles (including subdirectories).
# A Makefile that is not listed here will not be created by 'automake'!
AC_CONFIG_FILES([Makefile
                 extlib/Makefile
		 extlib/ape/Makefile
		 extlib/ape/include/Makefile
		 extlib/ape/include/ape/Makefile
		 extlib/ape/src/Makefile
		 extlib/atFunctions/Makefile
		 extlib/atFunctions/include/Makefile
		 extlib/atFunctions/src/Makefile
		 extlib/heaio/Makefile
		 extlib/heainit/Makefile
		 extlib/heautils/Makefile
		 extlib/heasp/Makefile
		 libsimput/Makefile
		 sample/Makefile
		 sample/csample/Makefile
		 sample/fsample/Makefile
		 tools/Makefile
		 tools/labnh/Makefile
		 tools/galabs/Makefile
		 tools/simputfile/Makefile
		 tools/simputlc/Makefile
		 tools/simputpsd/Makefile
		 tools/simputmerge/Makefile
		 tools/simputspec/Makefile
		 tools/simputsrc/Makefile])
		
AC_CONFIG_SUBDIRS([extlib/cfitsio extlib/wcslib extlib/fftw])

AC_OUTPUT

