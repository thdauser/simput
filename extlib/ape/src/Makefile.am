AM_CFLAGS=-I@top_srcdir@/extlib/ape/include -I@top_srcdir@/extlib/cfitsio -DUSE_READLINE

########### LIBRARIES #################

# Libraries that are installed in the 'lib' directory.
lib_LTLIBRARIES=libape.la

# Sources:
libape_la_SOURCES=ape_binary.c ape_error.c ape_io.c ape_list.c \
                  ape_msg.c ape_par.c ape_par_group.c ape_session.c ape_test.c \
                  ape_trad.c ape_util.c pil.c pil_error.c

libape_la_LDFLAGS = -version-info 1:0:0
