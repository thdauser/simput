#ifndef COMMON_H
#define COMMON_H 1


/////////////////////////////////////////////////////////////////
// Includes.
/////////////////////////////////////////////////////////////////


#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <assert.h>
#include <ctype.h>

#include "simput.h"


/////////////////////////////////////////////////////////////////
// Macros.
/////////////////////////////////////////////////////////////////


/** Common string length. */
#define SIMPUT_MAXSTR (1025)


/** Chatter level:

    0: error messages only

    1: error messages and warnings

    2: error messages, warnings, and informational output. 

    The chatter level can be defined as a compiler flag. If it is not
    set externally, use the default value of 1. */
#ifndef DCHATTY
#define DCHATTY 1
#endif

/** Output routine for error messages. */
#define SIMPUT_ERROR(msg) \
  fprintf(stderr, "\nError in %s: %s!\n", __func__, msg)

/** Chatter routine for warnings printed to STDOUT. */
#define SIMPUT_WARNING(msg) \
  if(1<=DCHATTY) {printf("\n*** Warning in %s: %s! ***\n", __func__, msg);}

/** Chatter routine for informational output to STDOUT. */
#define SIMPUT_INFO(msg) \
  if(2==DCHATTY) {printf("%s\n", msg);}


#define CHECK_STATUS_RET(a,b) \
  if (EXIT_SUCCESS!=a) return(b)

#define CHECK_STATUS_VOID(a) \
  if (EXIT_SUCCESS!=a) return

#define CHECK_STATUS_BREAK(a) \
  if (EXIT_SUCCESS!=a) break;

#define CHECK_NULL_RET(a,status,msg,ret) \
  if (NULL==a) { \
    SIMPUT_ERROR(msg); \
    status=EXIT_FAILURE; \
    return(ret);\
  }

#define CHECK_NULL_VOID(a,status,msg) \
  if (NULL==a) { \
    SIMPUT_ERROR(msg); \
    status=EXIT_FAILURE; \
    return;\
  }

#define CHECK_NULL_BREAK(a,status,msg) \
  if (NULL==a) { \
    SIMPUT_ERROR(msg); \
    status=EXIT_FAILURE; \
    break;\
  }

/** Macro returning the maximum of 2 values. */
#define MAX(a, b) ( (a)>(b) ? (a) : (b) )

/** Macro returning the minimum of 2 values. */
#define MIN(a, b) ( (a)<(b) ? (a) : (b) )

/** The following macros are used to the store light curve and the PSD
    in the right format for the GSL routines. */
#define REAL(z,i) ((z)[(i)])
#define IMAG(z,i,n) ((z)[(n)-(i)])


#endif /* COMMON_H */
