/*
 * userio.c: i/o routines for error, warning, log and progress messages
 */
#include <time.h>

#include "R.h"

#include "defs.h"

#include "debug.h"
#include "utils.h"
#include "s.h"
#include "userio.h"

static const char *error_messages[MAX_ERRNO+1] = {
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
/* 13 */	"out of dynamic memory (try local kriging?)",
/* 14 */	"i/o error: %s",
/* 15 */	"no command file%s",
/* 16 */	"%s user interface not compiled in this version",
/* 17 */	"writing to pipe `%s' failed",
/* 18 */	"reading from pipe `%s' failed",
/* 19 */    "function call prevented by secure mode%s"
};

/*
 * error handling function -- print message and error to string, and
 * call error message handler.
 */
void gstat_error(char *fname, int line, enum Gstat_errno err_nr, const char *msg) {

	assert(err_nr <= MAX_ERRNO);
	if (DEBUG_DUMP || err_nr == ER_NULL) /* print file&line */
		Rprintf("(%s, line %d)", fname, line);

	if (err_nr == ER_NULL)
		Rf_error("NULL error: this indicates a bug, please consider reporting this\n");

	if (msg == NULL)
		Rf_error("<NULL> message: indicating a software bug, please report\n");
	else
		Rf_error(error_messages[err_nr], msg);
	return;
}

/* message() calls for messages preceding a call to ErrMsg() */
void message(char *fmt, ...) {
	va_list args;
	/* char *buf = NULL; */
	char w[ERROR_BUFFER_SIZE];

	w[0] = '\0';
	va_start(args, fmt);
	vsnprintf(w, ERROR_BUFFER_SIZE, fmt, args);
	va_end(args);
	Rprintf("%s", w);
}

/* print a warning message to string, and call warning message handler */
void pr_warning(char *fmt, ...) {
	va_list args;
	char w[ERROR_BUFFER_SIZE];
	if (DEBUG_SILENT)
		return;

	w[0] = '\0';
	va_start(args, fmt);
	vsnprintf(w, ERROR_BUFFER_SIZE, fmt, args);
	va_end(args);
	Rf_warning("%s\n", w);
}

void printlog(const char *fmt, ...) {
	va_list args;
	char w[ERROR_BUFFER_SIZE];
	if (DEBUG_SILENT)
		return;

	w[0] = '\0';
	va_start(args, fmt);
	vsnprintf(w, ERROR_BUFFER_SIZE, fmt, args);
	va_end(args);
	Rprintf("%s", w);
}

void print_progress(unsigned int current, unsigned int total) {
	static int perc_last = -1, sec_last = -1;
	int perc, sec;
	static time_t start;

	R_CheckUserInterrupt(); /* allow for user interrupt */
	if (total <= 0 || DEBUG_SILENT || ! do_print_progress)
		return;

	if (sec_last == -1) {
		start = time(NULL);
		sec_last = 0;
	}
	perc = floor(100.0 * current / total);
	if (perc != perc_last) { /* another percentage -> calculate time: */
		if (current == total) { /* 100% done, reset: */
			Rprintf("\r%3d%% done\n", 100);
			perc_last = sec_last = -1;
		} else {
			sec = difftime(time(NULL), start);
			if (sec != sec_last) { /* another second -- don't print too often */
				Rprintf("\r%3d%% done", perc);
				perc_last = perc;
				sec_last = sec;
			}
		}
	}
}
