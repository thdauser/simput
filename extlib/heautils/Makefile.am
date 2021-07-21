AM_CFLAGS =-I@top_srcdir@/extlib/cfitsio
AM_CFLAGS+=-I@top_srcdir@/extlib/ape/include

########### LIBRARIES #################

# Libraries that are installed in the 'lib' directory.
lib_LTLIBRARIES=libhdutils.la

# Sources:
libhdutils_la_SOURCES=headas_copykeys.c headas_error.c \
                      headas_error_init.c headas_error_manager.c headas_error_map.c \
                      headas_history.c headas_toolname.c headas_utils.c \
                      headas_parstamp.c headas_gti.c \
                      HDgtcalf.c HDgtcalf_internal.c HDgtcalf_f77_wrap.c \
                      headas_file_check.c headas_rand.c mt.c headas_lock.c \
                      headas_polyfit.c headas_smooth.c headas_sort.c \
                      headas_svdfit.c

libhdutils_la_LDFLAGS = -version-info 1:0:0

libhdutils_la_LIBADD=../ape/src/libape.la

############ HEADERS #################

# Define the headers that will be installed in $(includedir):
include_HEADERS=headas_copykeys.h headas_error.h \
                headas_utils.h headas_gti.h hdcal.h headas_rand.h headas_lock.h \
                hd_math_compat.h headas_polyfit.h
