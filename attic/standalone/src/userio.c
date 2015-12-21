/*
    Gstat, a program for geostatistical modelling, prediction and simulation
    Copyright 1992, 2011 (C) Edzer Pebesma

    Edzer Pebesma, edzer.pebesma@uni-muenster.de
	Institute for Geoinformatics (ifgi), University of Münster 
	Weseler Straße 253, 48151 Münster, Germany. Phone: +49 251 
	8333081, Fax: +49 251 8339763  http://ifgi.uni-muenster.de 

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version. As a special exception, linking 
    this program with the Qt library is permitted.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

    (read also the files COPYING and Copyright)
*/

/*
 * userio.c: i/o routines for error, warning, log and progress messages
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <time.h>
#include <math.h>

#include "defs.h"

#ifdef PCRCALC
# define efclose fclose
#else
# include "err.h"
#endif

#include "debug.h"
#include "utils.h"
#include "version.h"
#include "userio.h"

#ifdef USING_R
void Rprintf(const char *, ...);
void Rf_error(const char *, ...);
# define is_openf(f) (f != NULL)
#else
# define is_openf(f) (f != NULL && f != stdout && f != stderr)
#endif


static FILE *logfile = NULL;

static struct {
	void (*warning_handler)(const char *mess);
	void (*error_handler)(const char *mess, int level);
	void (*printlog_handler)(const char *mess);
	void (*progress_handler)(unsigned int this, unsigned int total);
} gstat_handler = { NULL, NULL, NULL, NULL };

static void (*old_progress_handler)(unsigned int this, unsigned int total) 
		= NULL;

static STRING_BUFFER 
	*error_prefix = NULL, 
	*error_message = NULL, 
	*warning_message = NULL;

static enum Gstat_errno gstat_errno;

const char *error_messages[MAX_ERRNO+1] = {
/* 0 */		"%s",
/* 1 */		"bug in function `%s'",
/* 2 */		"variable not set: %s",
/* 3 */		"variable outside valid range: %s",
/* 4 */		"value not allowed for: %s",
/* 5 */		"no filename set %s",
/* 6 */		"write failed on file `%s'",
/* 7 */		"read failed on file `%s'",
/* 9 */		"cannot read real value from `%s'",
/* 9 */		"cannot read integer from `%s'",
/* 10 */	"syntax error: %s",
/* 11 */	"illegal option or missing argument on `%s'",
/* 12 */	"domain (math) error on `%s'",
/* 13 */	"out of dynamic memory %s",
/* 14 */	"i/o error: %s",
/* 15 */	"no command file%s",
/* 16 */	"%s user interface not compiled in this version",
/* 17 */	"writing to pipe `%s' failed",
/* 18 */	"reading from pipe `%s' failed",
/* 19 */    "function call prevented by secure mode%s",
/* 20 */	"matrix library error: %s",
/* 21 */	"extdbase error: %s"
};

void init_userio(int use_stdio) {
	if (use_stdio) {
#ifdef USING_R
		set_gstat_log_file(NULL);
#else
		set_gstat_log_file(stdout);
#endif
		set_gstat_warning_handler(default_warning);
		set_gstat_error_handler(default_error);
		set_gstat_log_handler(default_printlog);
		set_gstat_progress_handler(default_progress);
	} else {
		/* ... */
	}
	error_prefix    = resize_strbuf(error_prefix, ERROR_BUFFER_SIZE);
	error_message   = resize_strbuf(error_message, ERROR_BUFFER_SIZE);
	warning_message = resize_strbuf(warning_message, ERROR_BUFFER_SIZE);
	error_prefix->str[0] = error_message->str[0] = 
			warning_message->str[0] = '\0';
}

/*
 * error handling function -- print message and error to string, and
 * call error message handler.
 */
void gstat_error(char *fname, int line,
	enum Gstat_errno err_nr, const char *msg) {
	char s[30], *buf;
	int len;

	assert(err_nr <= MAX_ERRNO);
	gstat_errno = err_nr;

	if (error_prefix->str[0] != '\0')
		save_strcat(error_message, error_prefix->str);

	save_strcat(error_message, "gstat: ");
	len = strlen(error_message->str);
	buf = error_message->str + len;
#ifdef HAVE_SNPRINTF
	snprintf(buf, ERROR_BUFFER_SIZE - len,
		error_messages[err_nr], save_string(msg));
#else
	sprintf(buf, error_messages[err_nr], save_string(msg));
#endif
	if (DEBUG_DUMP || err_nr == ER_NULL) { /* print file&line */
		save_strcat(error_message, " (");
		save_strcat(error_message, fname);
		sprintf(s, ", line %d)", line);
		save_strcat(error_message, s);
	}

	if (err_nr == ER_NULL) {
		save_strcat(error_message, "\nVersion info: ");
		save_strcat(error_message, GSTAT_OS);
		save_strcat(error_message, " ");
		save_strcat(error_message, VERSION);
		save_strcat(error_message,
			"\nThis is a bug. Please send the above information, along with\n");
		save_strcat(error_message,
			"the information necessary to reproduce this bug to ");
		save_strcat(error_message, GSTAT_EMAIL);
	}

	gstat_handler.error_handler(error_message->str, err_nr);
	error_message->str[0] = '\0';
	return;
}

/* wrapper function for ErrClo(optopt), in case of error command line option */
void gstat_clo_error(char *f, int l, enum Gstat_errno err, int a) {
	static char s[2];
	sprintf(s, "%c", a);
	gstat_error(f, l, err, s);
} 

/* message() calls for messages preceding a call to ErrMsg() */
void message(char *fmt, ...) {
	va_list args;
	/* char *buf = NULL; */

	va_start(args, fmt);
#ifdef HAVE_VSNPRINTF
	vsnprintf(error_prefix->str, ERROR_BUFFER_SIZE, fmt, args);
#else
	vsprintf(error_prefix->str, fmt, args);
#endif
	va_end(args);
	/* buf = NULL; */
}

/* print a warning message to string, and call warning message handler */
void pr_warning(char *fmt, ...) {
	va_list args;
	char *buf = NULL;

	if (warning_message->max_length < 11)
		resize_strbuf(warning_message, 11);

	warning_message->str[0] = '\0';
	save_strcat(warning_message, "Warning: ");

	buf = warning_message->str + 9;

	va_start(args, fmt);
#ifdef HAVE_VSNPRINTF
	vsnprintf(buf, ERROR_BUFFER_SIZE - 9, fmt, args);
#else
	vsprintf(buf, fmt, args);
#endif
	va_end(args);

	gstat_handler.warning_handler(warning_message->str);
}

void print_progress(unsigned int current, unsigned int total) {
	gstat_handler.progress_handler(current, total);
}

/* get the value of gstat errno */
enum Gstat_errno get_gstat_errno(void) {
	return gstat_errno;
}

/* set the internal gstat errno to NO_ERROR, and reset error mesages */
void reset_gstat_errno(void) {
	assert(error_prefix);
	assert(error_message);

	gstat_errno = ER_NOERROR;
	error_prefix->str[0] = '\0';
	error_message->str[0] = '\0';
}

void set_gstat_warning_handler(void (*warning_fn)(const char *message)) {
	gstat_handler.warning_handler = warning_fn;
}

void set_gstat_error_handler(void (*error_fn)(const char *message, int level)) {
	gstat_handler.error_handler = error_fn;
}

void set_gstat_log_handler(void (*logprint)(const char *str)) {
	gstat_handler.printlog_handler = logprint;
}

void set_gstat_progress_handler(
		void (*progress)(unsigned int this, unsigned int total)) {
	gstat_handler.progress_handler = progress;
}

void push_gstat_progress_handler(
		void (*progress)(unsigned int this, unsigned int total)) {

	assert(old_progress_handler == NULL);

	old_progress_handler = gstat_handler.progress_handler;
	set_gstat_progress_handler(progress);
}

void pop_gstat_progress_handler(void) {

	assert(old_progress_handler != NULL);

	set_gstat_progress_handler(old_progress_handler);
	old_progress_handler = NULL;
}

const char *get_gstat_error_message(void) {
	return (const char *) error_message->str;
}

void print_to_logfile_if_open(const char *mess) {

	if (is_openf(logfile))
#ifdef USING_R
		Rprintf("%s", mess);
#else
		fprintf(logfile, "%s", mess);
#endif 
}

void default_warning(const char *mess) {

	print_to_logfile_if_open(mess);

#ifdef USING_R
	Rprintf("%s\n", mess);
#else
	fprintf(stderr, "%s\n", mess);
#endif
	return;
}

void default_error(const char *mess, int level) {

	print_to_logfile_if_open(mess);

#ifdef USING_R
	Rf_error("%s\n", mess);
#else
	fprintf(stderr, "%s\n", mess);
	exit(level == 0 ? -1 : level);
#endif
}

void printlog(const char *fmt, ...) {
	STRING_BUFFER *s;
	va_list args;

	s = resize_strbuf(NULL, ERROR_BUFFER_SIZE);

	va_start(args, fmt);
#ifdef HAVE_VSNPRINTF
	vsnprintf(s->str, ERROR_BUFFER_SIZE, fmt, args);
#else
	vsprintf(s->str, fmt, args);
#endif
	va_end(args);

	gstat_handler.printlog_handler(s->str);
	free_strbuf(s);
}

void default_printlog(const char *mess) {

	if (DEBUG_SILENT)
		return;

	if (is_openf(logfile))
		print_to_logfile_if_open(mess);
	else
#ifndef USING_R
		Rprintf("%s", mess);
#else
		printf("%s", mess);
#endif
}

int set_gstat_log_file(FILE *f) {
	int retval;

	if (f == NULL) {
		if (is_openf(logfile)) {
			retval = efclose(logfile);
			logfile = NULL;
			return retval;
		} else {
			logfile = NULL;
			return 1;
		}
	} else
		logfile = f;
	return 0;
}

void default_progress(unsigned int current, unsigned int total) {
	static int perc_last = -1, sec_last = -1;
	int perc, sec;
	static time_t start;

	if (total <= 0 || DEBUG_SILENT)
		return;

	if (sec_last == -1) {
		start = time(NULL);
		sec_last = 0;
	}
	perc = floor(100.0 * current / total);
	if (perc != perc_last) { /* another percentage -> calculate time: */
		if (current == total) { /* 100% done, reset: */
#ifdef USING_R
			Rprintf("\r%3d%% done\n", 100);
#else
			fprintf(stderr, "\r%3d%% done\n", 100);
#endif
			perc_last = sec_last = -1;
		} else {
			sec = difftime(time(NULL), start);
			if (sec != sec_last) { /* another second -- don't print too often */
#ifdef USING_R
				Rprintf("\r%3d%% done", perc);
#else
				fprintf(stderr, "\r%3d%% done", perc);
#endif
				perc_last = perc;
				sec_last = sec;
			}
		}
	}
}

#ifndef PCRCALC
/**************************** meschach error functions ****************/

#define SING_ERR \
"Read the manual at http://www.gstat.org/ ;\n\
look for: Trouble shooting -> Error messages -> From meschach"

#define MEM_ERR \
"In case you are trying to do global kriging (i.e., no neighbourhood\n\
parameters like `radius' or `max' were specified) with a large data set,\n\
reduce the neighbourhood size and use local kriging. In case you are\n\
fitting a variogram model to a large data set with REML, try another\n\
fitting method."

#define FP_ERR \
"This error may arise from using _very_ large or very small values for\n\
data values, variograms or coordinates. Try to rescale them to a\n\
reasonable range."

void setup_meschach_error_handler(int using_R) {
 	int code;
 	char *err, *hint = "", buf[100];

#ifndef PCRCALC
 	/* set up meschach error handler: */
 	if ((code = setjmp(restart)) == 0) {
		set_err_flag(using_R ? EF_R_ERROR /* avoid longjmp */
			: EF_JUMP /* make meschach jump on errors */  );
 	} else {
 		/* setjmp() returned non-zero, so we returned from a longjmp(): */
 		switch (code) {
 			case E_MEM:  /* run out of memory */
 				err = "virtual memory exhausted";
 				hint = MEM_ERR;
 				break;
 			case E_SING:  /* singular matrix occurred */
 				err = "singular matrix";
 				hint = SING_ERR;
 				break;
 			case E_POSDEF: /* non-positive definite matrix */
 				err = "non-positive definite matrix";
 				hint = "";
 				break;
 			case E_SIGNAL: /* floating point signal */
 				err = "floating point exception";
 				hint = FP_ERR;
 				break;
 			default:
 				sprintf(buf, "unknown, error code %d", code);
 				err = buf;
 				hint = "";
 				break;
 		}
		printlog("\ngstat caught an error that occurred in the matrix library,\n");
 		printlog("the reason for it was: %s\n\n", err);
		if (*hint)
			printlog("HINT: %s\n\n", hint);
 		ErrMsg(ER_MESCHACH, err);
 	} 
#endif /* PCRCALC */
}
#endif

void no_progress(unsigned int current, unsigned int total) {
	return;
}
