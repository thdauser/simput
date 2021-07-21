#include <stdio.h>
#include <stdlib.h>
#include <pil.h>
#include "hdcal.h"
#include "HDgtcalf_internal.h"

#ifdef HDGTCALF_STANDALONE
#include "HDgtcalf_standalone.h"
#else
#include "headas_error.h"
#endif

/*
int HDgtcalf (const char* tele, const char* instr, const char* detnam, 
      const char* filt, const char* codenam, const char* strtdate,
      const char* strtime, const char* stpdate, const char* stptime,
      const char* expr, int maxret,int fnamesize, char** filenam,
      long* extno, char** online, int* nret, int* nfound, int* status);
*/
int main() {

  int status = HD_OK;
  int func_status = 137;

  long extno = 0;
  int nret = 0;
  int nfound = 0;
  int maxret = 1;

  char *filename[1];
  char *online[1];

  filename[0]=(char *)malloc((sizeof(char))*PIL_LINESIZE);
  online[0]=(char *)malloc((sizeof(char))*32);

  func_status = HDgtcalf(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &func_status);
  if (137 != func_status) {
    fprintf(stderr, "HDgtcalf called with a non-0 status returned %d, not 137\n", func_status);
    status = 1;
  }

  func_status = 0;

  HDgtcalf(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &func_status);
  if (HD_ERR_NULL_POINTER != func_status) {
    fprintf(stderr, "HDgtcalf called with some 0s returned %d, not HD_ERR_NULL_POINTER (%d)\n", func_status,
      HD_ERR_NULL_POINTER);
    status = 1;
  }

  func_status = HDgtcalf(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
  if (HD_ERR_NULL_POINTER != func_status) {
    fprintf(stderr, "HDgtcalf called with all 0s returned %d, not HD_ERR_NULL_POINTER (%d)\n", func_status,
      HD_ERR_NULL_POINTER);
    status = 1;
  }

  /* Test sample input (Swift clock correction file).  Test may fail if CALDB 
     variables are not set correctly. */
  func_status = 0;

  HDgtcalf("swift", "SC", "-", "-", "CLOCK", "now", "now", "now", "now", "", maxret, 1024, filename, &extno, online, &nret, &nfound, &func_status);
  if (0 != func_status) {
    fprintf(stderr, "HDgtcalf called with sample input returned %d, not 0\n", func_status);
    /*
    fprintf(stderr, "HDgtcalf returned filename=%s\n", filename[0]);
    fprintf(stderr, "HDgtcalf returned nfound=%d\n", nfound);
    fprintf(stderr, "HDgtcalf returned nret=%d\n", nret);
    */
    status = 1;
  }

  return status;
}
