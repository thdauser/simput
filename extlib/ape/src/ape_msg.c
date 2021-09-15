/** \file ape_msg.c
    \brief Implementation of message facilities.
    \author James Peachey, HEASARC/EUD/GSFC.
*/
#include "ape/ape_msg.h"
#include "ape/ape_test.h"

#include <stdarg.h>
#include <stdio.h>
#include <string.h>

#define MSG_BUF_SIZE (FILENAME_MAX + 128)

#ifdef __cplusplus
extern "C" {
#endif

/* Stream to which error messages are directed. */
static FILE * s_stderr = 0;

/* Stream to which output is directed. */
static FILE * s_stdout = 0;

/* Flag indicating whether debug mode is active. */
static char s_debug_mode = 0;

/* Get error message stream, by default stderr, or else some other stream the client is using. */
static FILE * error_stream(void) {
  if (0 == s_stderr) s_stderr = stderr;
  return s_stderr;
}

static FILE * out_stream(void) {
  if (0 == s_stdout) s_stdout = stdout;
  return s_stdout;
}

static void write_err(const char * msg) {
  fprintf(error_stream(), "%s", msg);
  fflush(error_stream());
}

static void write_out(const char * msg) {
  fprintf(out_stream(), "%s", msg);
  fflush(out_stream());
}

static void write_warn(int chatter, const char * msg) {
  fprintf(error_stream(), "%s", msg);
  fflush(error_stream());
}

static void write_info(int chatter, const char * msg) {
  fprintf(out_stream(), "%s", msg);
  fflush(out_stream());
}

static ApeMsgFuncPtrType s_err_handler = &write_err;

static ApeMsgFuncPtrType s_out_handler = &write_out;

static ApeMsgChatPtrType s_warn_handler = &write_warn;

static ApeMsgChatPtrType s_info_handler = &write_info;

/* Forward a debugging message to error stream. */
void ape_msg_debug(const char * fmt, ...) {
  if (0 != s_debug_mode) {
    char msg[MSG_BUF_SIZE] = "";
    va_list ap;
    va_start(ap, fmt);
    vsprintf(msg, fmt, ap);
    va_end(ap);
    (*s_err_handler)(msg);
    va_end(ap);
  }
}

void ape_msg_debug_enable(char enable) { s_debug_mode = enable; }

char ape_msg_get_debug_mode(void) { return s_debug_mode; }

/* Forward an error message to error stream. */
void ape_msg_error(const char * fmt, ...) {
  char msg[MSG_BUF_SIZE] = "";
  va_list ap;
  va_start(ap, fmt);
  vsprintf(msg, fmt, ap);
  va_end(ap);
  (*s_err_handler)(msg);
}

/* Forward an informational message to error stream. */
/* TODO implement standard chatter facility. */
void ape_msg_info(unsigned int chatter, const char * fmt, ...) {
  char msg[MSG_BUF_SIZE] = "";
  va_list ap;
  va_start(ap, fmt);
  vsprintf(msg, fmt, ap);
  va_end(ap);
  (*s_info_handler)(chatter, msg);
}

void ape_msg_out(const char * fmt, ...) {
  char msg[MSG_BUF_SIZE] = "";
  va_list ap;
  va_start(ap, fmt);
  vsprintf(msg, fmt, ap);
  va_end(ap);
  (*s_out_handler)(msg);
}

/* Forward a warning message to error stream. */
/* TODO implement standard chatter facility. */
void ape_msg_warn(unsigned int chatter, const char * fmt, ...) {
  char msg[MSG_BUF_SIZE] = "";
  va_list ap;
  va_start(ap, fmt);
  vsprintf(msg, fmt, ap);
  va_end(ap);
  (*s_warn_handler)(chatter, msg);
}

/* Return the current ape error stream. */
FILE * ape_msg_get_err_stream(void) { return error_stream(); }

/* Reset the error stream to point to the new_stream. If new_stream is 0, de facto error messages will be
   redirected to stderr.  */
void ape_msg_set_err_stream(FILE * new_stream) {
  s_stderr = new_stream;
  /* Attempt to make stream unbuffered. */
  setbuf(s_stderr, 0);
}

ApeMsgFuncPtrType ape_msg_get_err_handler(void) { return s_err_handler; }

void ape_msg_set_err_handler(ApeMsgFuncPtrType func) {
  if (0 == func) s_err_handler = &write_err;
  else s_err_handler = func;
}

ApeMsgChatPtrType ape_msg_get_warn_handler(void) { return s_warn_handler; }

void ape_msg_set_warn_handler(ApeMsgChatPtrType func) {
  if (0 == func) s_warn_handler = &write_warn;
  else s_warn_handler = func;
}

/* Return the current ape output stream. */
FILE * ape_msg_get_out_stream(void) { return out_stream(); }

/* Reset the output stream to point to the new_stream. If new_stream is 0, de facto output messages will be
   redirected to stdout.  */
void ape_msg_set_out_stream(FILE * new_stream) { s_stdout = new_stream; }

ApeMsgFuncPtrType ape_msg_get_out_handler(void) { return s_out_handler; }

void ape_msg_set_out_handler(ApeMsgFuncPtrType func) {
  if (0 == func) s_out_handler = &write_out;
  else s_out_handler = func;
}

ApeMsgChatPtrType ape_msg_get_info_handler(void) { return s_info_handler; }

void ape_msg_set_info_handler(ApeMsgChatPtrType func) {
  if (0 == func) s_info_handler = &write_info;
  else s_info_handler = func;
}

/* Display a unit test failure. Note this is really part of the ape_test interface but it appears here so
   that it can mimic ape_msg_error. */
void ape_test_failed(const char * fmt, ...) {
  char msg[MSG_BUF_SIZE] = "";
  va_list ap;
  va_start(ap, fmt);
  strcat(msg, "ERROR: ");
  vsprintf(msg + strlen(msg), fmt, ap);
  va_end(ap);
  (*s_err_handler)(msg);
  /* In addition to displaying the message, set the global unit test status to indicate a test failed. */
  ape_test_set_status(1);
}

/* A test-writing function which does nothing with its input. */
void null_writer(const char * msg) {
}

void null_chatter(int chatter, const char * msg) {
}

/* Unit test of message facilities. */
void ape_msg_test() {
  /* Save state of debugging flag. */
  char debug_mode = s_debug_mode;

  /* Turn off debugging mode for the next test. */
  ape_msg_debug_enable(0);

  /* Write a message to debugging stream before enabling debugging, using one argument. */
  ape_msg_debug("THIS DEBUG MESSAGE SHOULD NEVER APPEAR!.\n", 1);

  /* Restore debugging mode state. */
  ape_msg_debug_enable(debug_mode);

  /* Write a message to debugging stream, using one argument. */
  ape_msg_debug("This debug message should always appear %d time in the output file.\n", 1);

  /* Write a message to error stream, using one argument. */
  ape_msg_error("This error message should always appear %d time in the output file.\n", 1);

  /* Write a message to warning stream, using chatter 0 (always displayed). */
  ape_msg_warn(0, "This warning should always appear %s\n", "in the output file.");

  /* Set output stream handler to use the null handler. */
  ape_msg_set_out_handler(null_writer);

  /* Write message which shouldn't be seen to test the redirect. */
  ape_msg_out("THIS OUTPUT MESSAGE SHOULD BE REDIRECTED TO NOWHERE, AND SHOULD NEVER APPEAR!.\n");

  /* Set output stream handler to 0 (default handler). */
  ape_msg_set_out_handler(0);

  /* Set info stream handler to use the null handler. */
  ape_msg_set_info_handler(null_chatter);

  /* Write message which shouldn't be seen to test the redirect. */
  ape_msg_info(3, "THIS INFO MESSAGE SHOULD BE REDIRECTED TO NOWHERE, AND SHOULD NEVER APPEAR!.\n");

  /* Set info stream handler to 0 (default handler). */
  ape_msg_set_info_handler(0);

  /* Set error stream handler to use the null handler. */
  ape_msg_set_err_handler(null_writer);

  /* Write message which shouldn't be seen to test the redirect. */
  ape_msg_error("THIS ERROR MESSAGE SHOULD BE REDIRECTED TO NOWHERE, AND SHOULD NEVER APPEAR!.\n");

  /* Set error stream handler to 0 (default handler). */
  ape_msg_set_err_handler(0);

  /* Set warn stream handler to use the null handler. */
  ape_msg_set_warn_handler(null_chatter);

  /* Write message which shouldn't be seen to test the redirect. */
  ape_msg_warn(3, "THIS WARNING MESSAGE SHOULD BE REDIRECTED TO NOWHERE, AND SHOULD NEVER APPEAR!.\n");

  /* Set warn stream handler to 0 (default handler). */
  ape_msg_set_warn_handler(0);

  /* Write message which should be seen to test the restore. */
  ape_msg_out("This message should go to the default output stream.\n");

  /* Write message which should be seen to test the restore. */
  ape_msg_info(3, "This message should go to the default info stream.\n");

  /* Write message which should be seen to test the restore. */
  ape_msg_error("This message should go to the default error stream.\n");

  /* Write message which should be seen to test the restore. */
  ape_msg_warn(3, "This message should go to the default warning stream.\n");

}

#ifdef __cplusplus
}
#endif

/*
 * $Log$
 * Revision 1.9  2013/09/06 19:14:49  peachey
 * Add ape_msg_set_err_handler function, parallel to ape_msg_set_out_handler.
 * The new function allows the client code to supply a custom error handling
 * function.
 *
 * Revision 1.8  2006/06/16 01:18:48  peachey
 * Add ape_msg_get_debug_mode, for getting current debug mode.
 *
 * Revision 1.7  2006/06/03 02:11:38  peachey
 * Allow redirection of output stream using functions.
 *
 * Revision 1.6  2006/05/19 17:36:55  peachey
 * Unbuffer redirected error stream.
 *
 * Revision 1.5  2006/05/19 17:33:43  peachey
 * Flush error stream; needed in case it was redirected to an
 * buffered stream.
 *
 * Revision 1.4  2006/05/19 17:30:34  peachey
 * Add ape_msg_set_out_stream, for redirecting stdout the way ape already
 * redirects stderr.
 *
 * Revision 1.3  2006/05/18 13:58:42  peachey
 * Add ape_msg_out, for output.
 *
 * Revision 1.2  2006/04/12 17:56:34  peachey
 * Add ERROR banner to ape_test_failed, and add ape_msg_info function.
 *
 * Revision 1.1.1.1  2006/04/05 13:45:19  peachey
 * Initial import of All-purpose Parameter Environment (APE).
 *
*/
