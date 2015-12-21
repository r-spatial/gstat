#ifndef USERIO_H
#define USERIO_H

enum Gstat_errno {
	ER_NOERROR     =  0 /* no error */,
	ER_NULL        =  1 /* internal error: should not occur */,
	ER_VARNOTSET   =  2 /* a required variable was not set by the user */,
	ER_RANGE       =  3 /* range error (outside permitted values) */,
	ER_IMPOSVAL    =  4 /* a variable was set to an illegal value */,
	ER_WRITE       =  6 /* write error on file */,
	ER_READ        =  7 /* read error on file */,
	ER_RDFLT       =  8 /* error while converting a string to a float */,
	ER_RDINT       =  9 /* error while converting a string to an int */,
	ER_SYNTAX      = 10 /* syntax error */,
	ER_ARGOPT      = 11 /* error in command line option arguments */,
	ER_DOMAIN      = 12 /* math error */,
	ER_MEMORY      = 13 /* memory exhausted */,
	ER_IO          = 14 /* i/o conflict (e.g. redirection not permitted) */,
	ER_NOCMD       = 15 /* no command file specified */,
	ER_NOCURSES    = 16 /* no curses user interface compiled in */,
	ER_PWRITE      = 17 /* error while writing to a pipe */,
	ER_PREAD       = 18 /* error while reading from a pipe */,
	ER_SECURE      = 19 /* secure mode: operation not allowed */
};

#define MAX_ERRNO 19
void message(char *fmt, ...);  /* message() calls always preceed ErrMsg() */
#define ErrMsg(a,b) gstat_error(__FILE__,__LINE__,a,b)
void gstat_error(char *fname, int line, enum Gstat_errno err_nr, const char *msg);

void pr_warning(char *fmt, ...);
void printlog(const char *fmt, ...);
void print_progress(unsigned int current, unsigned int total);

#endif /* USERIO_H */
