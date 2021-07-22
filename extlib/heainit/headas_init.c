/******************************************************************************
 *   File name: headas_init.c                                                 *
 *                                                                            *
 * Description: Universal initialization code for any HEADAS task. Normally   *
 *     this is called by headas_main.c, and need not be called explicitly.    *
 *                                                                            *
 *    Language: C or C++                                                      *
 *                                                                            *
 *      Author: Mike Tripicco, for HEASARC/GSFC/NASA                          *
 *                                                                            *
 *  Change log: see CVS Change log at the end of the file.                    *
 ******************************************************************************/

/******************************************************************************
 * Header files.                                                              *
 ******************************************************************************/

#ifdef AHLOG_INTEGRATION_COMPLETE_REMOVE_THIS_IF
#include "ahlog/cahlog.h"
#endif

#include "ape/ape_error.h"
#include "ape/ape_msg.h"
#include "ape/ape_trad.h"

#include <ctype.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "headas.h"
#include "headas_error.h"
#include "pil.h"
/******************************************************************************/

/* C/C++ compatibility. */
#ifdef __cplusplus
extern "C" {
#endif

  /****************************************************************************
   * Static variable definitions.                                             *
   ****************************************************************************/
static char sPILInitSuccess = 0;
  /****************************************************************************/

  /****************************************************************************
   * Static function definitions.                                             *
   ****************************************************************************/
static int StandardPfile(int argc, char *argv[]);
static int hd_pil_err_logger(char *s);
static void ReportError(int taskStatus);
  /****************************************************************************/

  /****************************************************************************
   * Function definitions.                                                    *
   ****************************************************************************/
/* called by headas_main.c prior to tool subroutine */
int headas_init(int argc, char *argv[]){
    int status = 0;
    int tmpstatus = 0;

    headas_chatpar = -1;
    headas_clobpar = 0;

    /* Register defaults for toolname/version, used in error reporting.
       These functions cannot create an error condition. */
    set_toolname(hdbasename(argv[0]));
    set_toolversion("0.0");

    /* read env vars and setup output streams */
    tmpstatus = HDIO_init();
    if (tmpstatus) HDerror_throw(0, 0, 0, tmpstatus);

    status = status ? status : tmpstatus;

    /* Set up global error handler. This can/should be done even if an
       error occurred with the streams. */
    tmpstatus = HDerror_init(0);
    if (tmpstatus) HDerror_throw(0, 0, 0, tmpstatus);

    status = status ? status : tmpstatus;

    tmpstatus = StandardPfile(argc, argv);

    status = status ? status : tmpstatus;

    return status;
}

/* called by headas_main.c after tool subroutine */
int headas_close(int taskStatus) {
  int status = 0;

  /* Close PIL provided the taskStatus argument was not a PIL error. */
  if(sPILInitSuccess &&
      (PIL_ERR_BASE < taskStatus || PIL_ERR_MIN_IDX > taskStatus)) {
    if(PIL_OK != (status = PILClose(PIL_OK))) {
      HD_ERROR_THROW("PILClose failed", status);
    }
  }

  /* If the taskStatus argument is 0 (no error before this function
     was called) then set it to the local error status. */
  taskStatus = taskStatus ? taskStatus : status;

  /* Report a non-0 overall task error status. */
  if(taskStatus && !HDerror_dump_is_silent()) ReportError(taskStatus);

  return taskStatus;
}

int StandardPfile(int argc, char *argv[]) {
    int status;
    int headas_histpar = -1;

    if (PIL_OK != (status = PILInit(argc, argv))){
	HD_ERROR_THROW("PIL initialization failed", status);
	return status;
    }
    sPILInitSuccess = 1;

/* enable 'batch mode' via HEADASNOQUERY environment variable */
    if (HD_getenv("HEADASNOQUERY")) PILOverrideQueryMode(PIL_QUERY_OVERRIDE);

    status = PILSetLoggerFunction(hd_pil_err_logger);
    if(status != 0){
        HD_ERROR_THROW("PILSetLoggerFunction failed", status);
        return status;
    }

    /* Use heautils's function to check for file existence. */
    status = PILSetFileAccessFunction(&HDfile_check);
    if(status != 0){
        HD_ERROR_THROW("PILSetFileAccessFunction failed", status);
        return status;
    }

    /* unset temporarily to keep missing parameters     */
    /* (chatter/clobber/history) from throwing an error */
    if (PIL_OK != (status = PILSetLoggerFunction(NULL))) {
        HD_ERROR_THROW("PILSetLoggerFunction failed", status);
        return status;
    }

    if (PIL_OK != (status=PILGetInt("chatter", &headas_chatpar))) {
      /* OK not to have chatter parameter     */
      /* but we'll leave headas_chatpar at -1 */
	status = 0;
    }

    status=PILGetBool("clobber", &headas_clobpar);
    if(status != 0){ /* OK not to have clobber parameter */
	status = 0;
    }

    status=PILGetBool("history", &headas_histpar);
    if(status != 0){ /* OK not to have history parameter     */
	             /* but we'll leave headas_histpar at -1 */
	status = 0;
    }

    /* now reset to normal logging mode */
    status = PILSetLoggerFunction(hd_pil_err_logger);
    if(status != 0){
        HD_ERROR_THROW("PILSetLoggerFunction failed", status);
        return status;
    }

    /* Handle history parameter. */
    set_history(headas_histpar);

    return status;
}

static int hd_pil_err_logger(char *s)
{
   fprintf(heaerr, "%s\n", s);
   return(0);
}

static void ReportError(int taskStatus) {
  char task[73] = "";
  char vers[9] = "";
  const char *code;
  const char *pkgid;
  const char *text;
  char found;
  char doDump;
  char cfitsioMsg[FLEN_ERRMSG];

  /* First, dump the Cfitsio stack. */
  doDump = 1;
  while(HD_OK != fits_read_errmsg(cfitsioMsg)) {
    if(doDump) {
      fprintf(heaerr, "Dumping CFITSIO error stack:\n--------------------------------------------------\n");
      doDump = 0; /* Only write the header the first time. */
    }

    fprintf(heaerr, "%s\n", cfitsioMsg);
  }
  /* If doDump is not 1, the while loop above printed something, so
     print a footer. */
  if(0 == doDump) fprintf(heaerr, "--------------------------------------------------\nCFITSIO error stack dump complete.\n");

  /* If HEADAS error stack does not contain an error, throw one now, so
     that HDerror_dump can print the error message. */
  if(0 == HDerror_get()) HDerror_throw(NULL, NULL, 0, taskStatus);

  /* Dump the HEADAS message stack. */
  HDerror_dump(heaerr);

  /* Get error information. */
  HDerror_get_info(taskStatus, &found, &pkgid, &code, &text);

  /* If no info was available, get it from the system. */
  if(!found && 0 < taskStatus) {
    text = strerror(taskStatus);
    if(text) found = 1;
  }

  /* Print error information summary. */
  if(found && (pkgid || code || text)) {
     if(pkgid) fprintf(heaerr, "%s ", pkgid);
     fprintf(heaerr, "ERROR");
     if(code) fprintf(heaerr, " %s", code);
     if(text) fprintf(heaerr, ": %s", text);
     fprintf(heaerr, "\n");
  }

  /* Get information about the task. */
  get_toolname(task);
  get_toolversion(vers);

  /* Print information about the task. This is the final message
     written to heaerr. */
  fprintf(heaerr, "Task");
  if('\0' != *task) fprintf(heaerr, " %s", task);
  if('\0' != *vers) fprintf(heaerr, " %s", vers);
  fprintf(heaerr, " terminating with status %d\n", taskStatus);

  /* Make certain the returned status is in the range the shell can
     handle. */
  taskStatus = 201; /* Generic headas return value. */

  /* Give specific codes for known components. */
  if(pkgid) {
    if(strstr(pkgid, "CFITSIO")) taskStatus = 202;
    else if(strstr(pkgid, "PIL")) taskStatus = 203;;
  }
}
  /****************************************************************************/

#ifdef AHLOG_INTEGRATION_COMPLETE_REMOVE_THIS_IF

#define HDMSGLEN (80)

static void ape_out_redirect(const char* msg) {
  ahlog_out(__func__, "%s", msg);
}

static void ape_err_redirect(const char* msg) {
  ahlog_warn(LOW, __func__, "%s", msg);
}

int headas_start_up(int argc, char * argv[], const char * tooltag) {
  /* Set up logging. Ahlog uses parameters for this: task name,
       logfile, chatter, debug. Heainit uses environment variables.
     Ahapp and heainit provide local functions to redirect ape logging.
     Both read/use clobber, history.
TODO Move the following comment block somewhere:
     Order: goal is to set up logging ASAP, then other things as required.
  */
/*  std::string banner(__func__); banner += ": "; */
  int chatter = 0;
  char chatter_warn_msg[HDMSGLEN] = "";

  char * logfilepar = 0;
  const char * logfile = 0;
  char logfile_warn_msg[HDMSGLEN] = "";

  char c_debug = 0;
  char debug_warn_msg[HDMSGLEN] = "";

  char clobber = 0;

  char c_history = 0;

  int status = ape_trad_init(argc, argv);
  if (eOK != status) {
    /* TODO: report status code along with message. */
    fprintf(stderr, "%s: %s %d.\n", __func__, "ape_trad_init returned status", status);
/*    throw std::runtime_error(banner + "ape_trad_init returned non-0 status"); */
    return status;
  }

  /* Disable all prompting if HEADASNOQUERY environment variable is set. */
  if (0 != getenv("HEADASNOQUERY")) PILOverrideQueryMode(PIL_QUERY_OVERRIDE);

  set_history(1);

  /* attempt to read chatter from parameter file
     ... if chatter out of range, adjust to appropriate limit
     ... if other error, set chatter to max value */
/*  int chatter=0; */
  /* std::string chatter_warn_msg=""; */
  status = ape_trad_query_int("chatter", &chatter);
  if (eValueBelowMin == status) {
    chatter=MINCHAT;
    /* std::stringstream msg; */
    sprintf(chatter_warn_msg, "Chatter too small; resetting to minimum value of %d", MINCHAT);
    /* chatter_warn_msg=msg.str(); */
  } else if (eValueAboveMax == status) {
    chatter=MAXCHAT;
    /* std::stringstream msg; */
    sprintf(chatter_warn_msg, "Chatter too large; resetting to maximum value of %d", MAXCHAT);
    /* chatter_warn_msg=msg.str(); */
  } else if (eOK != status) {
    chatter=MAXCHAT;
    /* std::stringstream msg; */
    sprintf(chatter_warn_msg, "Invalid chatter value; resetting to maximum value of %d", MAXCHAT);
    /* chatter_warn_msg=msg.str(); */
  }
  status=0;    /* reset ape status */

  /* attempt to read log file name */
  /* ... if ape error, set log file name to !DEFAULT */
  /* char* logfilepar = 0;
  const char* logfile = 0;
  std::string logfile_warn_msg=""; */
  status = ape_trad_query_string("logfile",&logfilepar);
  if (eOK == status && strlen(logfilepar) > 0) {
    logfile = logfilepar;
  } else {
    logfile = "!DEFAULT"; /* Is this really what we want: create log file if logfile parameter missing? */
    strcat(logfile_warn_msg, "APE error reading logfile parameter; setting logfile to !DEFAULT");
  }
  status=0;    /* reset ape status */

  /* attempt to read debug parameter */
  /* ... if ape error, set debug parameter to false */
  /* char c_debug = 0;
  bool debug=false;
  std::string debug_warn_msg=""; */
  status = ape_trad_query_bool("debug", &c_debug);
  if (eOK != status) {
    /* debug=false; */
    strcat(debug_warn_msg, "APE error reading debug parameter; setting debug to false");
  } else {
    /* debug = 0 != c_debug ? true : false; */
  }
  status=0;    /* reset ape status */

  /* set up ahlog as soon as possible to have as many messages sent through ahlog */
  ahlog_setup(argv[0], logfile, chatter, c_debug);
  free(logfilepar); logfilepar=0;

  /* display any necessary warning messages */
  if ('\0' != *chatter_warn_msg) ahlog_warn(LOW, __func__, "%s\n", chatter_warn_msg);
  if ('\0' != *logfile_warn_msg) ahlog_warn(LOW, __func__, "%s\n", logfile_warn_msg);
  if ('\0' != *debug_warn_msg) ahlog_warn(LOW, __func__, "%s\n", debug_warn_msg);

  /* feed APE error messages through ahlog */
  ape_msg_set_out_handler(&ape_out_redirect);
  ape_msg_set_err_handler(&ape_err_redirect);

  /* attempt to read clobber parameter */
  /*  ... if ape error, set clobber to false */
  /* char clobber = 0; */
  /* bool b_clobber=false; */
  status = ape_trad_query_bool("clobber", &clobber);
  if (eOK != status) {
    headas_clobpar=0;
    ahlog_warn(LOW, __func__,  "APE error reading clobber parameter; setting clobber to false\n");
  } else {
    /* b_clobber = 0 != clobber ? true : false; */
    headas_clobpar = clobber;
  }
  status=0;
  /* ahgen::setClobber(b_clobber); */ /* TODO: need to find a suitable replacement? */
  ahlog_set_clobber(clobber);
#if 0
  /* Skip the following block forever -- it's ahfits-specific. */
  /* attempt to read buffer parameter */
  /*  ... if buffer < -1; set to buffer=-1 */
  /*  ... if ape error; set buffer=-1 */
  int buffer = -1;    /* default is automatic buffering */
  status = ape_trad_query_int("buffer", &buffer);
  if (eOK != status) {
    buffer=-1;
    AH_WARN(ahlog::LOW) << "APE error reading buffer parameter; setting buffer to -1 (automatic buffering)" << std::endl;
  } else if (buffer < -1) {
    buffer=-1;
    AH_WARN(ahlog::LOW) << "buffer parameter set smaller than -1; resetting buffer to -1 (automatic buffering)" << std::endl;
  }
  status=0;
  ahgen::setBuffer(buffer);
#endif

  /* attempt to read history parameter which enables/disables parameter stamping */
  /* char c_history = 0; */
  /* bool history=true; */
  status = ape_trad_query_bool("history", &c_history);
  if (eOK != status) {
    /* history=true; */
    set_history(1);
    ahlog_warn(LOW, __func__, "%s\n", "APE error reading history parameter; setting history to true");
  } else {
    /* history = 0 != c_history ? true : false; */
    set_history(c_history);
  }
  status=0;
  /* ahgen::setHistory(history); */ /* TODO: need to find a suitable replacement? */

  /* make sure that stamped parameters match those used; this is necessary since */
  /* in case ape errors above cause a hard-coded default value to be used */
  ape_trad_set_int("chatter",chatter);
  ape_trad_set_string("logfile",logfile);
  ape_trad_set_bool("debug",c_debug);
  ape_trad_set_bool("clobber",clobber);
  /* ape_trad_set_int("buffer",buffer); */
  ape_trad_set_bool("history",c_history);

  /* set tool name used when stamping parameters */
  /* std::string tname=argv[0]; */
  /* size_t found; */
  /* found=tname.rfind("/"); */
  /* if (found != std::string::npos) tname.replace(0,found+1,""); */
  { char * tname = strrchr(argv[0], '/');
    if (0 == tname) tname = argv[0];
    else ++tname; /* Go past the slash. */
    set_toolname(tname);
    ahlog_set_executable_name(tname);
  }

  /* set tool version */
  /* the tool version is typically set using a macro set as the CVS macro, */
  /* Name, so that the tooltag argument to this function will take the form: */
  /* $Name$.  Before setting the tool version, we want to strip */
  /* the CVS formatting, namely "$Name$" at the */
  /* end. */
  /* std::string ttooltag=tooltag; */
  /* std::string tagstr="$Name: "; */
  /* note: in the following search, the leading $ is ommitted to prevent CVS from expanding the keyword */
  /* if (ttooltag.find("Name$") != std::string::npos) ttooltag=""; */     /* check if no CVS flag set yet */
  /* if (ttooltag[ttooltag.size()-1] == '$') */
  /*   ttooltag.erase(ttooltag.end()-1,ttooltag.end()); */
  /* if (ttooltag.find(tagstr) != std::string::npos) */
  /*   ttooltag.replace(ttooltag.find(tagstr),tagstr.size(),""); */
  /* set_toolversion(ttooltag.c_str()); */
  { const char * begin = strstr(tooltag, "$Name");
    const char * end = 0;
    char * ttooltag = 0;
    size_t ttooltaglen = 0u;
    if (0 == begin) {
      begin = tooltag;
    } else {
      /* Skip over the initial $Name. */
      begin += strlen("$Name");

      /* Skip the ":", if it is present. */
      if (':' == *begin) ++begin;

      /* Remove any leading space. */
      while ('\0' != *begin && 0 != isspace(*begin)) ++begin;

      /* Terminate the string with the next occurrence of $. */
      end = strchr(begin, '$');

      /* Remove trailing space, if any. */
      while (end > begin && 0 != isspace(*(end-1))) --end;
    }
    /* Make sure end has been set, and find length of version tag. */
    if (0 == end) end = begin + strlen(begin);
    ttooltaglen = end - begin;

    /* Copy the significant portion of the tag into ttooltag string. */
    ttooltag = (char *) calloc(ttooltaglen + 1, sizeof(char)); /* Include room for terminator. */
    if (0 == ttooltag) status = 1;
    else strncpy(ttooltag, begin, ttooltaglen);

    /* Store the tool version string and clean up. */
    if (0 != ttooltag) set_toolversion(ttooltag);
    free(ttooltag); ttooltag = 0;
  }
  return status;
}

int headas_shut_down(void) {
  ape_trad_close(1);
  /* TODO: more shut-down needed here? */
  return 0;
}
#endif /* AHLOG_INTEGRATION_COMPLETE_REMOVE_THIS_IF */

/* C/C++ compatibility. */
#ifdef __cplusplus
}
#endif

/******************************************************************************
 * $Log$
 * Revision 1.19  2014/01/03 22:31:47  peachey
 * Add and test new functionality headas_start_up.
 *
 * Revision 1.18  2006/06/16 00:21:01  peachey
 * Suppress attempting to print system error messages if error code is negative.
 *
 * Revision 1.17  2006/06/07 07:09:26  peachey
 * Do not call deprecated function PILVerifyCmdLine; PILInit does this.
 *
 * Revision 1.16  2004/08/17 19:30:57  peachey
 * Use new HDfile_check function to check for existence of files.
 *
 * Revision 1.15  2004/01/07 14:46:37  peachey
 * Check whether errors have been silenced before reporting them.
 *
 * Revision 1.14  2003/02/13 20:03:50  peachey
 * Refine error reporting.
 *
 * Revision 1.13  2003/02/13 18:45:55  peachey
 * Perform both error and stream initializations even if an error
 * occurs. This allows meaningful error reporting in such cases.
 *
 * Revision 1.12  2003/02/13 12:25:02  peachey
 * Do not close headas streams; they will be closed automatically on exit.
 * Otherwise, spurious errors can crop up from streams being closed more
 * than once.
 *
 * Revision 1.11  2003/02/12 23:27:45  peachey
 * hdIOInit code moved to headas_stdio.c, and renamed HDIO_init.
 *
 * Revision 1.10  2002/11/06 19:15:26  peachey
 * Initialize heain stream.
 *
 * Revision 1.9  2002/10/29 17:12:57  peachey
 * Do not put PIL error messages on the HEADAS stack; instead write them
 * directly to heaerr stream.
 *
 * Revision 1.8  2002/10/18 19:39:25  peachey
 * Reduce the number of dashes on the CFITSIO error stack banner.
 *
 * Revision 1.7  2002/10/09 16:26:39  peachey
 * Improve error handling in init/close. Use streams heaout, heaerr, heaprom
 * throughout. Move some code out of higher level functions into lower
 * level functions.
 *
 * Revision 1.6  2002/10/04 21:51:03  peachey
 * Changes to use the new error handling facility automatically.
 *
 ******************************************************************************/
