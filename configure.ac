#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.65)
AC_INIT([simput],[0.33],[christian.schmid@sternwarte.uni-erlangen.de])
# Store auxiliary files in the sub-directory 'build-aux'.
AC_CONFIG_AUX_DIR([build-aux]) 
AM_INIT_AUTOMAKE([foreign -Wall -Werror])
AC_PROG_CC # Search for a C compiler and store it in the Variable 'CC'.
AC_PROG_FC # Search for a Fortran compiler.
AC_CONFIG_SRCDIR([libsimput/fileaccess.c])

LT_INIT #([disable-shared]) #([shared]) 

# Compiler and Linker Flags:
# "-g -W -Wall -O0" for the compiler for debugging
# "-g" to create debugging information (seems to be used automatically by autotools)
# "-Werror" for the compiler to treat warnings as errors
# "-pg" for the compiler and linker to include information for gprof (Profiler)
# "-std=c99 -D_GNU_SOURCE" for C99 syntax standards
# " -O0" for zero optimization
CFLAGS="$CFLAGS -g -std=c99 -D_GNU_SOURCE -W -Wall -fPIC"
LDFLAGS="$LDFLAGS -g -W -Wall"
# FCFLAGS="-O0"

####################################
# Conditional compilation of Fortran code.
AC_ARG_ENABLE([fortran], 
              [AS_HELP_STRING([--disable-fortran], 
	                          [disable compilation of Fortran code])],
	          [fortran=false], [fortran=true])
AM_CONDITIONAL([FORTRAN], [test x$fortran = xtrue])

####################################
# Check for libraries:  (objdump -R/-T libncurses...)
AC_SEARCH_LIBS(sin, m, [], [AC_MSG_ERROR([No math library found!])], [])

#################################### 
# Check for header files:
AC_HEADER_STDC
AC_CHECK_HEADERS([limit.h math.h stdio.h stdlib.h stdarg.h string.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST

# The following line tells 'configure' to create a header file 'config.h'
# that is filled with definitions from this script.
# The 'config.h' contains for example the name of the package and can be
# included in the programm sources.
AC_CONFIG_HEADER([config.h])

# Create the following Makefiles (including subdirectories).
# A Makefile that is not listed here will not be created by 'automake'!
AC_CONFIG_FILES([Makefile
                 extlib/Makefile
		 libsimput/Makefile
		 sample/Makefile
		 sample/csample/Makefile
		 sample/fsample/Makefile])
		
AC_CONFIG_SUBDIRS([extlib/cfitsio extlib/wcslib extlib/fftw])

AC_OUTPUT

